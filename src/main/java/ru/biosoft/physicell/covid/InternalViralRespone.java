package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Phenotype;

public class InternalViralRespone
{

    //    void internal_virus_response_model_setup(  )
    //    {
    //        // set up the model 
    //            // set version info 
    //        internal_virus_response_model_info.name = "internal viral response"; 
    //        internal_virus_response_model_info.version = internal_virus_response_version; 
    //            // set functions 
    //        internal_virus_response_model_info.main_function = NULL; 
    //        internal_virus_response_model_info.phenotype_function = internal_virus_response_model; 
    //        internal_virus_response_model_info.mechanics_function = NULL; 
    //        
    //            // what pCell.getMicroenvironment() variables do you expect? 
    //            // what custom data do I need? 
    //        internal_virus_response_model_info.cell_variables.push_back( "max_infected_apoptosis_rate" ); 
    //        internal_virus_response_model_info.cell_variables.push_back( "max_apoptosis_half_max" ); 
    //        internal_virus_response_model_info.cell_variables.push_back( "apoptosis_hill_power" );  
    //        
    //            // register the submodel  
    //        internal_virus_response_model_info.register_model();    
    //            // set functions for the corresponding cell definition 
    //            
    //    //  pCD = pCell.getModel().getCellDefinition( "lung epithelium" ); 
    //    //  pCD.functions.update_phenotype = epithelium_submodel_info.phenotype_function;
    //    //  pCD.functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
    //        
    //        return; 
    //    }

    static void internal_virus_response_model(Cell pCell, Phenotype phenotype, double dt)
    {
        CellDefinition pCD = pCell.getModel().getCellDefinition( "lung epithelium" );

        // bookkeeping -- find pCell.getMicroenvironment() variables we need

        int nV_external = pCell.getMicroenvironment().findDensityIndex( "virion" );
        int nA_external = pCell.getMicroenvironment().findDensityIndex( "assembled virion" );
        int chemokine_index = pCell.getMicroenvironment().findDensityIndex( "chemokine" );
        int nINF1 = pCell.getMicroenvironment().findDensityIndex( "interferon 1" );

        int nV_internal = pCell.customData.findVariableIndex( "virion" );
        int nA_internal = pCell.customData.findVariableIndex( "assembled_virion" );
        int nP = pCell.customData.findVariableIndex( "viral_protein" );


        // actual model goes here 

        // now, set apoptosis rate 

        int apoptosis_model_index = pCell.phenotype.death.findDeathModelIndex( "apoptosis" );
        // phenotype.death.rates[apoptosis_model_index] = 

        // base death rate (from cell line)
        double base_death_rate = pCD.phenotype.death.rates.get( apoptosis_model_index );

        // additional death rate from infectoin  
        double additional_death_rate = pCell.customData.get( "max_infected_apoptosis_rate" );


        double v = pCell.customData.get( nA_internal ) / pCell.customData.get( "max_apoptosis_half_max" );
        v = Math.pow( v, pCell.customData.get( "apoptosis_hill_power" ) );

        double effect = v / ( 1.0 + v );
        additional_death_rate *= effect;
        phenotype.death.rates.set( apoptosis_model_index, base_death_rate + additional_death_rate );

        // if we're infected, secrete a chemokine for the immune model
        double AV = pCell.customData.get( nA_internal );

        /* old 
        if( P > 0.001 )
        {
            phenotype.secretion.secretionRates[chemokine_index] = 
                pCell.customData.get( "infected_cell_chemokine_secretion_rate" ];
            phenotype.secretion.saturation_densities[chemokine_index] = 1.0;        
        }
        else
        {
            phenotype.secretion.secretionRates[chemokine_index] = 0.0;
        }
        */

        // if I am dead, make sure to still secrete the chemokine 

        //  int chemokine_index = pCell.getMicroenvironment().findDensityIndex( "chemokine" ); 
        //  int nP = pCell.customData.findVariableIndex( "viral_protein"); 
        // double P = pCell.customData.get(nP];
        int proinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "pro-inflammatory cytokine" );

        int nR = pCell.customData.findVariableIndex( "viral_RNA" );
        double R = pCell.customData.get( nR );
        int antibody_index = pCell.getMicroenvironment().findDensityIndex( "Ig" );



        if( R >= 1.00 - 1e-16 )
        {
            pCell.customData.set( "infected_cell_chemokine_secretion_activated", 1.0 );
            // (AJ-V5) Antibody binding starts once the cell is infected
            pCell.phenotype.secretion.uptakeRates[antibody_index] = pCell.getModel().getParameterDouble( "Antibody_binding_rate" );
        }

        if( pCell.customData.get( "infected_cell_chemokine_secretion_activated" ) > 0.1 && phenotype.death.dead == false )
        {
            double rate = AV;
            rate /= pCell.customData.get( "max_apoptosis_half_max" );
            if( rate > 1.0 )
            {
                rate = 1.0;
            }
            rate *= pCell.customData.get( "infected_cell_chemokine_secretion_rate" );

            phenotype.secretion.secretionRates[chemokine_index] = rate;
            phenotype.secretion.saturationDensities[chemokine_index] = 1.0;

            // (Adrianne) adding pro-inflammatory cytokine secretion by infected cells
            pCell.phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = pCell.customData
                    .get( "activated_cytokine_secretion_rate" );
        }

        // (Adrianne) check whether the cell is undergoing pyroptosis and if so, evaluate the pyropotosis model
        if( pCell.customData.get( "cell_pyroptosis_flag" ) > 0.5 )
        {
            pyroptosis_cascade( pCell, phenotype, dt );
            return;
        }
        ;

        // Sara&Fiona: Internalise pro-pyroptotic cytokine from the pCell.getMicroenvironment()
        int pro_pyroptotic_cytokine = pCell.getMicroenvironment().findDensityIndex( "pro-pyroptosis cytokine" );
        double pyroptotic_cytokine_concentration = phenotype.molecular.internSubstrates[pro_pyroptotic_cytokine];

        //printf("PCC: %lf\n",pyroptotic_cytokine_concentration);

        //phenotype.secretion.uptake_rates[propyroptotic_cytokine_index] = 1.; // What should the secretion intake depend on ?
        //double pyroptotic_cytokine_concentration = phenotype.molecular.internalized_total_substrates[propyroptotic_cytokine_index];
        //phenotype.secretion.uptake_rates[nV_external] = pCell.customData.get(nR_bind] * pCell.customData.get(nR_EU]; 

        // (Sara&Fiona) pyroptosis cascade in the cell is initiated if cell's viral_RNA is >1 (i.e. >=3). This is arbitraty to check things work.
        if( R >= 150 && (int)pCell.customData.get( "cell_pyroptosis_flag" ) == 0
                && (int)pCell.customData.get( "cell_virus_induced_apoptosis_flag" ) == 0 )
        {
            // set the probability (in 0,1) that a cell with a death-sentence pyroptoses (not apoptoses)
            double cell_death_pyroptosis_probability = ( R - 100 ) / ( 1000 - 100 );
            if( cell_death_pyroptosis_probability > 1.0 )
            {
                cell_death_pyroptosis_probability = 1.0;
            }
            cell_death_pyroptosis_probability /= 2;
            // randomise a number in 0,1 that determines the cell death mode (pyroptosis or apoptosis)
            if( pCell.getModel().getRNG().UniformRandom() < cell_death_pyroptosis_probability )
            {
                pCell.customData.set( "cell_pyroptosis_flag", 1 ); //cell pyroptoses
            }

            return;
        }
        // (Sara&Fiona)
        else if( pyroptotic_cytokine_concentration > 100.0 && (int)pCell.customData.get( "cell_pyroptosis_flag" ) == 0 )
        {
            pCell.customData.set( "cell_pyroptosis_flag", 1 ); // Pyroptosis cascade is initiated
            pCell.customData.set( "cell_bystander_pyroptosis_flag", 1 ); // Pyroptosis cascade is initiated
            //printf("Pyro bystander effect!\n");
            return;
        }

        // interferon signaling 

        // // approximate activation by extracellular interferons 
        // // // activation = min( 1 , extracellular_interferon / interferon_max_response_threshold ) 
        pCell.customData.set( "interferon_activation",
                pCell.nearest_density_vector()[nINF1] / ( pCell.customData.get( "interferon_max_response_threshold" ) + 1e-32 ) );
        if( pCell.customData.get( "interferon_activation" ) > 1.0 )
        {
            pCell.customData.set( "interferon_activation", 1.0 );
        }

        // // Type-I interferon secretion 
        // // // secretion_rate = r_viral * Heaviside( RNA - 1 ) + r_paracrine * activation 
        phenotype.secretion.secretionRates[nINF1] = pCell.customData.get( "interferon_activation" );
        phenotype.secretion.secretionRates[nINF1] *= pCell.customData.get( "max_interferon_secretion_rate_via_paracrine" );
        if( R >= pCell.customData.get( "interferon_viral_RNA_detection" ) - 1e-16 ) // if there is at least 1 complete set of uncoated viral RNA
        {
            double scaled_RNA = ( R - pCell.customData.get( "interferon_viral_RNA_detection" ) )
                    / ( -pCell.customData.get( "interferon_viral_RNA_detection" )
                            + pCell.customData.get( "interferon_viral_RNA_threshold" ) );
            if( scaled_RNA > 1 )
            {
                scaled_RNA = 1.0;
            }
            phenotype.secretion.secretionRates[nINF1] += pCell.customData.get( "interferon_secretion_rate_via_infection" ) * scaled_RNA;
        }

        // // now the interferon response 
        // // // protein_synthesis_rate = protein_synthesis_rate_0 * ( 1 - interferon_activation * interferon_max_virus_inhibition ) 
        double protSynthRate = pCell.customData.get( "interferon_max_virus_inhibition" );
        protSynthRate *= pCell.customData.get( "interferon_activation" );
        protSynthRate *= -1;
        protSynthRate += 1.0;
        protSynthRate *= pCD.custom_data.get( "protein_synthesis_rate" );
        pCell.customData.set( "protein_synthesis_rate", protSynthRate ); // inhibition
        //        pCell.customData.set("protein_synthesis_rate") *= pCell.customData.get("interferon_activation"); // activation*inhibition
        //        pCell.customData.set("protein_synthesis_rate") *= -1; // -activation*inhibition 
        //        pCell.customData.set("protein_synthesis_rate") += 1.0; // 1 - activation*inhibition 
        //        pCell.customData.set("protein_synthesis_rate") *= pCD.customData.get("protein_synthesis_rate"); 
        // protein_synthesis_rate0 * (1 - activation*inhibition)  
    }

    // Sara&Fiona: code for pyroptosis
    static void pyroptosis_cascade(Cell pCell, Phenotype phenotype, double dt)
    {
        // Sara&Fiona: Pyroptosis code starts here.

        // Intracellular components
        int nfkb_n = pCell.customData.findVariableIndex( "nuclear_NFkB" );
        int nlrp3_i = pCell.customData.findVariableIndex( "inactive_NLRP3" );
        int nlrp3_a = pCell.customData.findVariableIndex( "active_NLRP3" );
        int nlrp3_b = pCell.customData.findVariableIndex( "bound_NLRP3" );
        int asc_b = pCell.customData.findVariableIndex( "bound_ASC" );
        int caspase1_b = pCell.customData.findVariableIndex( "bound_caspase1" );
        int gsdmd_c = pCell.customData.findVariableIndex( "cleaved_gasderminD" );
        int il_1b_p = pCell.customData.findVariableIndex( "pro_IL_1b" );
        int il_1b_c = pCell.customData.findVariableIndex( "cytoplasmic_IL_1b" );
        int il_1b_e = pCell.customData.findVariableIndex( "external_IL_1b" );
        int il_18_c = pCell.customData.findVariableIndex( "cytoplasmic_IL_18" );
        int il_18_e = pCell.customData.findVariableIndex( "external_IL_18" );
        int volume_c = pCell.customData.findVariableIndex( "cytoplasmic_volume" );
        //Not needed components (can be implicitely found via conservation laws)
        // double nfkb_c = pCell.customData.findVariableIndex( "cytoplasmic_NFkB_fraction" ); 
        // double asc_f = pCell.customData.findVariableIndex( "free_ASC" );
        // double caspase1_f = pCell.customData.findVariableIndex( "free_caspase1" );
        // double gsdmd_uc = pCell.customData.findVariableIndex( "uncleaved_gasderminD" );
        // double il_18_e = pCell.customData.findVariableIndex( "external_IL_1b" );

        //We can add these to the xml later..
        /**
        // Rate constants
        double k_nfkb_ctn = pCell.customData.findVariableIndex( "rate_NFkB_cytoplasm_to_nucleus" ); 
        double k_nfkb_ntc = pCell.customData.findVariableIndex( "rate_NFkB_nucleus_to_cytoplasm" );
        double k_nlrp3_ita = pCell.customData.findVariableIndex( "rate_NLRP3_incactive_to_active" );
        double k_nlrp3_atb = pCell.customData.findVariableIndex( "rate_NLRP3_active_to_bound" );
        double k_asc_ftb = pCell.customData.findVariableIndex( "rate_ASC_free_to_bound" );
        double k_c1_ftb = pCell.customData.findVariableIndex( "rate_caspase1_free_to_bound" );
        double k_il1b_cte = pCell.customData.findVariableIndex( "rate_Il1b_cytoplasmic_to_external" );
        double k_il18_cte = 1;//pCell.customData.findVariableIndex( "rate_Il18_cytoplasmic_to_external" );
        double k_vol_c = 1;// pCell.customData.findVariableIndex( "rate_pyroptosis_volume_increase" );
        // Decay constants
        double d_nlrp3 = pCell.customData.findVariableIndex( "decay_NLRPR_inactive_and_active" );
        double d_il = pCell.customData.findVariableIndex( "decay_IL1b" );
        // Hill function rates
        double a_nlrp3 = pCell.customData.findVariableIndex( "rate_constant_NLRP3_production" );
        double a_il1b_p = pCell.customData.findVariableIndex( "rate_constant_IL1b_production" );
        double a_gsdmd = pCell.customData.findVariableIndex( "rate_constant_GSDMD_cleavage" );
        double a_il1b_c = pCell.customData.findVariableIndex( "rate_constant_Il1b_cleavage" );
        double a_il18 = pCell.customData.findVariableIndex( "rate_constant_Il18_cleavage" );
        double hm_nfkb = pCell.customData.findVariableIndex( "halfmax_NFkB_transcription" );
        
        //std::cout<<pCell.customData.findVariableIndex( "halfmax_caspase1_cleavage" )<<std::endl;
        //double hm_c1 = 1;//pCell.customData.findVariableIndex( "halfmax_caspase1_cleavage" );
        double hm_c1 = pCell.customData.findVariableIndex("halfmax_caspase1_cleavage");
        double hex_nfkb = pCell.customData.findVariableIndex( "hillexponent_NFkB_transcription" );
        double hex_c1 = pCell.customData.findVariableIndex( "hillexponent_caspase1_cleavage" );
        // Total concentrations (let's have them all as 1 now and compute fractions)
        // double tot_nfkb = pCell.customData.findVariableIndex( "total_NFkB_concentration" );    
        // double tot_asc = pCell.customData.findVariableIndex( "total_ASC_concentration" );
        // double tot_c1 = pCell.customData.findVariableIndex( "total_caspase1_concentration" );
        // double tot_gsdmd = pCell.customData.findVariableIndex( "total_GSDMD_concentration" );
        // double tot_il18 = pCell.customData.findVariableIndex( "total_IL18_concentration" );
        
        */


        // System of pyroptosis equations starts here
        //Model constants (definitions could be moved to xml file)
        double k_nfkb_ctn = 0.3;
        double k_nfkb_ntc = 0.03;
        double k_nlrp3_ita = 0.07;
        double k_nlrp3_atb = 0.07;
        double k_asc_ftb = 0.02;
        double k_c1_ftb = 0.04;
        double k_il1b_cte = 0.8;
        double k_il18_cte = 0.8;
        double k_vol_c = 0.1;
        // Decay constants
        double d_nlrp3 = 0.002;
        double d_il = 0.004;
        // Hill function rates
        double a_nlrp3 = 0.025;
        double a_il1b_p = 0.007;
        double a_gsdmd = 0.08;
        double a_il1b_c = 0.8;
        double a_il18 = 0.8;
        double hm_nfkb = 0.3;
        double hm_c1 = 0.3;
        double hex_nfkb = 2.0;
        double hex_c1 = 2.0;
        //If the inflammsome base is formed set F_ib = 0. 
        double F_ib = 1;
        if( pCell.customData.get( nlrp3_b ) >= 1 )
        {
            F_ib = 0;
        }

        //Update nuclear NFkB (updated backward)
        pCell.customData.set( nfkb_n,
                ( pCell.customData.get( nfkb_n ) + k_nfkb_ctn * F_ib * dt ) / ( 1 + dt * k_nfkb_ntc + k_nfkb_ctn * F_ib * dt ) );

        //Set Hill function 1
        double hill_nfkb = ( Math.pow( pCell.customData.get( nfkb_n ), hex_nfkb ) )
                / ( Math.pow( hm_nfkb, hex_nfkb ) + Math.pow( pCell.customData.get( nfkb_n ), hex_nfkb ) );

        //Update NLRP3 (inactive, active and bound) (updated backward)
        pCell.customData.set( nlrp3_i,
                ( pCell.customData.get( nlrp3_i ) + dt * a_nlrp3 * hill_nfkb ) / ( 1 + dt * k_nlrp3_ita + dt * d_nlrp3 ) );
        pCell.customData.set( nlrp3_a, ( pCell.customData.get( nlrp3_a ) + k_nlrp3_ita * dt * ( pCell.customData.get( nlrp3_i ) ) )
                / ( 1 + dt * k_nlrp3_atb + dt * d_nlrp3 ) );
        pCell.customData.set( nlrp3_b, pCell.customData.get( nlrp3_b ) + dt * k_nlrp3_atb * F_ib * pCell.customData.get( nlrp3_a ) );

        //Update bound ASC (updated backward)
        pCell.customData.set( asc_b, ( pCell.customData.get( asc_b ) + dt * k_asc_ftb * ( 1 - F_ib ) * ( pCell.customData.get( nlrp3_b ) ) )
                / ( 1 + dt * k_asc_ftb * ( 1 - F_ib ) * ( pCell.customData.get( nlrp3_b ) ) ) );

        //Update bound caspase1 (updated backward)
        pCell.customData.set( caspase1_b, ( pCell.customData.get( caspase1_b ) + dt * k_c1_ftb * ( pCell.customData.get( asc_b ) ) )
                / ( 1 + dt * k_c1_ftb * ( pCell.customData.get( asc_b ) ) ) );

        //Set Hill function 2
        double hill_caspase1 = ( Math.pow( pCell.customData.get( caspase1_b ), hex_c1 ) )
                / ( Math.pow( hm_c1, hex_c1 ) + Math.pow( pCell.customData.get( caspase1_b ), hex_c1 ) );

        //Update cleaved GSDMD (updated backward)
        pCell.customData.set( gsdmd_c,
                ( pCell.customData.get( gsdmd_c ) + dt * a_gsdmd * hill_caspase1 ) / ( 1 + dt * a_gsdmd * hill_caspase1 ) );

        //Set G function (same now that total GSDMD concentration is 1 au of concentration)
        double g_gsdmd = pCell.customData.get( gsdmd_c ) / 1;

        //Update IL1b (pro, cytoplasmic, external)  We want to relate this to secreted cytokine IL1b (updated backward)
        pCell.customData.set( il_1b_p,
                ( pCell.customData.get( il_1b_p ) + dt * a_il1b_p * hill_nfkb ) / ( 1 + dt * a_il1b_c * hill_caspase1 + dt * d_il ) );
        pCell.customData.set( il_1b_c,
                ( pCell.customData.get( il_1b_c ) + dt * a_il1b_c * hill_caspase1 * ( pCell.customData.get( il_1b_p ) ) )
                        / ( 1 + dt * d_il + dt * k_il1b_cte * g_gsdmd ) );
        pCell.customData.set( il_1b_e, pCell.customData.get( il_1b_e ) + dt * ( k_il1b_cte * g_gsdmd * pCell.customData.get( il_1b_c ) ) );

        //Update IL18 (cytoplasmic, external)(updated backward)
        pCell.customData.set( il_18_c,
                ( pCell.customData.get( il_18_c ) + dt * a_il18 * hill_caspase1 * ( 1 - pCell.customData.get( il_18_e ) ) )
                        / ( ( 1 + dt * a_il18 * hill_caspase1 ) * ( 1 + dt * k_il18_cte * g_gsdmd ) ) );
        pCell.customData.set( il_18_e, pCell.customData.get( il_18_e ) + dt * k_il18_cte * g_gsdmd * pCell.customData.get( il_18_c ) );

        //Update cytoplasmic volume (updated backward)
        pCell.customData.set( volume_c, pCell.customData.get( volume_c ) / ( 1 - dt * k_vol_c * g_gsdmd ) );

        // (Yafei) need to update the real radius 
        phenotype.volume.total = pCell.customData.get( volume_c );

        //Temporary: "super fast" apoptosis occurs when cell should burst. 
        //To do: We actually want the cell to rupture once a cytoplasmic critical volume is reached (e.g. 1.5 of initial cytoplasmic volume from in vitro data). 
        int apoptosis_model_index = pCell.phenotype.death.findDeathModelIndex( "apoptosis" );
        double initial_total_volume = 2494;

        if( pCell.customData.get( volume_c ) > 1.2 * initial_total_volume )
        {
            //std::cout<<"Pyroptotic cell burst!"<<std::endl;
            //The cell's 'apoptosis death rate' is set to be "super high" 
            phenotype.death.rates.set( apoptosis_model_index, 9e9 );
        }
        // (Adrianne) update cell pro-inflammatory secretion rate based on IL18 secretion rate - need to double check unit conversion
        int proinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "pro-inflammatory cytokine" );
        pCell.phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = pCell.customData
                .get( "activated_cytokine_secretion_rate" ) + k_il18_cte * g_gsdmd * pCell.customData.get( il_18_c );
        // (Sara and Fiona)
        int propyroptotic_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "pro-pyroptosis cytokine" );
        pCell.phenotype.secretion.secretionRates[propyroptotic_cytokine_index] = k_il1b_cte * g_gsdmd * pCell.customData.get( il_1b_c );
    }

}
