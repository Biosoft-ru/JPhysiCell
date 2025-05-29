package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class epithelium_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {

        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        // receptor dynamics 
        // requires faster time scale - done in main function 

        // viral dynamics model 
        InternalVirusDynamics.internal_virus_model( pCell, phenotype, dt );
        //            internal_viral_dynamics_info.phenotype_function( pCell, phenotype, dt );
        // internal_virus_model(pCell,phenotype,dt);

        // viral response model 
        InternalViralRespone.internal_virus_response_model( pCell, phenotype, dt );
        //            internal_virus_response_model_info.phenotype_function( pCell, phenotype, dt );
        // internal_virus_response_model(pCell,phenotype,dt);   

        // T-cell based death
        TCell_induced_apoptosis( pCell, phenotype, dt );

        // (Adrianne V5) ROS induced cell death model
        ROS_induced_apoptosis( pCell, phenotype, dt );

        // if I am dead, remove all adhesions 
        int apoptosis_index = phenotype.death.findDeathModelIndex( "apoptosis" );
        if( phenotype.death.dead )
        {
            // detach all attached cells 
            // remove_all_adhesions( pCell ); 

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
        }

        /*
        // cell secretion belongs in viral response 
        
        // if I am dead, make sure to still secrete the chemokine 
        static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
        static int nP = pCell.custom_data.find_variable_index( "viral_protein"); 
        double P = pCell.custom_data[nP];
        
        static int nAV = pCell.custom_data.find_variable_index( "assembled_virion" ); 
        double AV = pCell.custom_data[nAV]; 
        
        static int nR = pCell.custom_data.find_variable_index( "viral_RNA");
        double R = pCell.custom_data[nR];
        
        if( R >= 1.00 - 1e-16 ) 
        {
        pCell.custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
        }
        
        if( pCell.custom_data["infected_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
        {
        double rate = AV; // P; 
        rate /= pCell.custom_data["max_apoptosis_half_max"];
        if( rate > 1.0 )
        { rate = 1.0; }
        rate *= pCell.custom_data[ "infected_cell_chemokine_secretion_rate" ];
        
        phenotype.secretion.secretion_rates[chemokine_index] = rate; 
        phenotype.secretion.saturation_densities[chemokine_index] = 1.0; 
        }
        */

        // if I am dead, don't bother executing this function again 
        if( phenotype.death.dead )
        {
            pCell.functions.updatePhenotype = null;
        }
    }

    void TCell_induced_apoptosis(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );
        int proinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "pro-inflammatory cytokine" );
        int antiinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "anti-inflammatory cytokine" );

        if( pCell.customData.get( "TCell_contact_time" ) > pCell.customData.get( "TCell_contact_death_threshold" ) )
        {
            // make sure to get rid of all adhesions! 
            // detach all attached cells 
            // remove_all_adhesions( pCell ); 

            //            #pragma omp critical
            //            {
            //            std::cout << "\t\t\t\t" << pCell << " (of type " << pCell.type_name <<  ") died from T cell contact" << std::endl; 
            //            }

            // induce death 
            pCell.startDeath( apoptosis_index );

            pCell.phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = 0;
            // pCell.phenotype.secretion.secretion_rates[antiinflammatory_cytokine_index] = pCell.custom_data["antiinflammatory_cytokine_secretion_rate_by_damagedSite"]; 
            pCell.phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
            ((ModelCovid)pCell.getModel()).addResidual(  pCell.position.clone() );
            pCell.functions.updatePhenotype = null;
        }
    }

    void ROS_induced_apoptosis(Cell pCell, Phenotype phenotype, double dt)
    {
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );
        int ROS_index = pCell.getMicroenvironment().findDensityIndex( "ROS" );
        double ROS_amount = pCell.nearest_density_vector()[ROS_index];
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );
        int proinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "pro-inflammatory cytokine" );

        double epsilon_ROS = pCell.getModel().getParameterDouble( "epsilon_ROS" );

        double prob_apoptosis = ROS_amount / ( ROS_amount + epsilon_ROS );

        if( pCell.getModel().getRNG().UniformRandom() < prob_apoptosis )
        {
            //            std::cout<<ROS_amount<<" "<<epsilon_ROS<<std::endl;
            // make sure to get rid of all adhesions! 
            // detach all attached cells 
            // remove_all_adhesions( pCell ); 

            //            #pragma omp critical
            //            {
            //            std::cout << "\t\t\t\t" << pCell << " (of type " << pCell.type_name <<  ") died from ROS" << std::endl; 
            //            }

            // induce death 
            pCell.startDeath( apoptosis_index );

            pCell.phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = 0;
            pCell.phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );

            pCell.functions.updatePhenotype = null;
        }
    }
}