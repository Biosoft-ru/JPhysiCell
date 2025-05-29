package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;

public class InternalVirusDynamics
{
    //    #include "./internal_viral_dynamics.h" 

    //    using namespace PhysiCell; 

    //    std::string internal_virus_replication_version = "0.5.0"; 

    //    Submodel_Information internal_viral_dynamics_info; 


    //    void internal_virus_model_setup(  )
    //    {
    //            // set version
    //        internal_viral_dynamics_info.name = "internal viral replication dynamics"; 
    //        internal_viral_dynamics_info.version = internal_virus_replication_version; 
    //            // set functions 
    //        internal_viral_dynamics_info.main_function = NULL; 
    //        internal_viral_dynamics_info.phenotype_function = internal_virus_model; 
    //        internal_viral_dynamics_info.mechanics_function = NULL; 
    //            // what microenvironment variables do I need? 
    //
    //            // what custom data do I need? 
    //        internal_viral_dynamics_info.microenvironment_variables.push_back( "assembled virion" );    
    //
    //        internal_viral_dynamics_info.cell_variables.push_back( "virion" ); // adhered, in process of endocytosis 
    //        internal_viral_dynamics_info.cell_variables.push_back( "uncoated_virion" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "viral_RNA" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "viral_protein" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "assembled_virion" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "export_virion" ); 
    //
    //        internal_viral_dynamics_info.cell_variables.push_back( "virion_uncoating_rate" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "uncoated_to_RNA_rate" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "protein_synthesis_rate" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "virion_assembly_rate" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "virion_export_rate" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "max_RNA_replication_rate" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "RNA_replication_half" ); 
    //        internal_viral_dynamics_info.cell_variables.push_back( "basal_RNA_degradation_rate" ); 
    //
    //        // submodel_registry.register_model( internal_viral_dynamics_info ); 
    //        internal_viral_dynamics_info.register_model();
    //        
    //        return; 
    //    }

    static void internal_virus_model(Cell pCell, Phenotype phenotype, double dt)
    {
        // bookkeeping -- find microenvironment variables we need

        int nV_external = pCell.getMicroenvironment().findDensityIndex( "virion" );
        int nA_external = pCell.getMicroenvironment().findDensityIndex( "assembled virion" );

        int nV_internal = pCell.customData.findVariableIndex( "virion" );
        int nA_internal = pCell.customData.findVariableIndex( "assembled_virion" );

        int nUV = pCell.customData.findVariableIndex( "uncoated_virion" );
        int nR = pCell.customData.findVariableIndex( "viral_RNA" );
        int nP = pCell.customData.findVariableIndex( "viral_protein" );
        int eP = pCell.customData.findVariableIndex( "export_virion" );

        /*  
         bool done = false; 
        extern Cell* pInfected; 
        if( pCell == pInfected && 1 == 0 )
        {
            std::cout << std::endl << "viral dynamics : " << __LINE__ << " " 
                << phenotype.molecular.internalized_total_substrates[ nV_external ] << " " 
                << phenotype.molecular.internalized_total_substrates[ nA_external ] << " " 
                << pCell.customData.get(nV_internal] << " " 
                << pCell.customData.get(nUV] << " " 
                << pCell.customData.get(nR] << " " 
                << pCell.customData.get(nP] << " "    
                << pCell.customData.get(nA_internal] << " " 
                << std::endl;       
        }
        */

        // copy virions from "internalized variables" to "custom variables"
        /*  
        pCell.customData.get(nV_internal] = 
            phenotype.molecular.internalized_total_substrates[nV_external]; 
        // this transfer is now handled in receptor dynamics 
        */

        // This resets the internal assembled virion count 
        // so we are commenting it out 
        // pCell.customData.get(nA_internal] = 
        //  phenotype.molecular.internalized_total_substrates[nA_external]; 

        // actual model goes here 

        // uncoat endocytosed virus
        double dV = dt * pCell.customData.get( "virion_uncoating_rate" ) * pCell.customData.get( nV_internal );
        if( dV > pCell.customData.get( nV_internal ) )
        {
            dV = pCell.customData.get( nV_internal );
        }
        pCell.customData.modify( nV_internal, -dV );
        pCell.customData.modify( nUV, dV );

        // convert uncoated virus to usable mRNA 
        double dR = dt * pCell.customData.get( "uncoated_to_RNA_rate" ) * pCell.customData.get( nUV );
        // gotta remove this from uncoated virions now befoe we add the replication
        if( dR > pCell.customData.get( nUV ) )
        {
            dR = pCell.customData.get( nUV );
        }
        pCell.customData.modify( nUV, -dR );
        // RNA replication post uncoated to RNA calc
        dR += dt * pCell.customData.get( "max_RNA_replication_rate" ) * pCell.customData.get( nR )
                / ( pCell.customData.get( nR ) + pCell.customData.get( "RNA_replication_half" ) );
        // RNA degradation
        dR -= dt * pCell.customData.get( "basal_RNA_degradation_rate" ) * pCell.customData.get( nR );

        if( dR < -1 * pCell.customData.get( nR ) )
        {
            dR = -1 * pCell.customData.get( nR );
        }
        pCell.customData.modify( nR, dR );

        // use mRNA to create viral protein 
        double dP = dt * pCell.customData.get( "protein_synthesis_rate" ) * pCell.customData.get( nR );
        pCell.customData.modify( nP, dP );

        // degrade protein 

        // assemble virus 
        double dA = dt * pCell.customData.get( "virion_assembly_rate" ) * pCell.customData.get( nP );
        if( dA > pCell.customData.get( nP ) )
        {
            dA = pCell.customData.get( nP );
        }
        pCell.customData.modify( nP, -dA );
        pCell.customData.modify( nA_internal, dA );

        // set export rate 
        /*  
        phenotype.secretion.net_export_rates[nA_external] = 
            pCell.customData.get("virion_export_rate" ] * pCell.customData.get(nA_internal]; 
         
        // copy data from custom variables to "internalized variables" 
            
        phenotype.molecular.internalized_total_substrates[nV_external] = 
            pCell.customData.get(nV_internal];        
        phenotype.molecular.internalized_total_substrates[nA_external] = 
            pCell.customData.get(nA_internal];    
        */
        double deP = dt * pCell.customData.get( "virion_export_rate" ) * pCell.customData.get( nA_internal );
        if( deP > pCell.customData.get( nA_internal ) )
        {
            deP = pCell.customData.get( nA_internal );
        }
        pCell.customData.modify( eP, deP );
        pCell.customData.modify( nA_internal, -deP );
        //test int export

        double alpha1 = Math.floor( pCell.customData.get( eP ) );

        pCell.nearest_density_vector()[nV_external] += alpha1 / pCell.getMicroenvironment().mesh.dV;
        pCell.customData.modify( eP, -alpha1 );
    }
}