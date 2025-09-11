package ru.biosoft.physicell.covid;

import java.util.Set;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class macrophage_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );
        CellDefinition pCD = pCell.getModel().getCellDefinition( "macrophage" );
        int proinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "pro-inflammatory cytokine" );
        int chemokine_index = pCell.getMicroenvironment().findDensityIndex( "chemokine" );
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );
        int virus_index = pCell.getMicroenvironment().findDensityIndex( "virion" );
        int antibody_index = pCell.getMicroenvironment().findDensityIndex( "Ig" );
        int antiinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "anti-inflammatory cytokine" );

        // no apoptosis until activation (resident macrophages in constant number for homeostasis) 
        if( pCell.customData.get( "activated_immune_cell" ) < 0.5 )
        {
            phenotype.death.rates.set( apoptosis_index, 0.0 );
        }
        else
        {
            phenotype.death.rates.set( apoptosis_index, pCD.phenotype.death.rates.get( apoptosis_index ) );
        }

        if( phenotype.death.dead == true )
        {
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
            return;
        }

        // make changes to volume change rate??

        // if too much debris, comit to apoptosis   

        /* // remove in v 3.2   
        double relative_volume = ( phenotype.volume.total/pCD.phenotype.volume.total ); 
        if( relative_volume > pCell.customData.get( "relative_maximum_volume" ] )
        {
            pCell.start_death( apoptosis_index ); 
            pCell.phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = 0; 
            pCell.phenotype.secretion.secretionRates[debris_index] = pCell.customData.get("debris_secretion_rate"]; 
            
            return;
        }
        */

        // check for cells to eat 
        Set<Cell> neighbors = pCell.cells_in_my_container();//  .cells_in_my_container(); 

        // at least one of the cells is pCell 
        if( neighbors.size() < 2 )
        {
            return;
        }

        // (Adrianne) get type of CD8+ T cell and CD4+ t CELL
        int CD8_Tcell_type = pCell.getModel().getCellDefinition( "CD8 Tcell" ).type;
        int CD4_Tcell_type = pCell.getModel().getCellDefinition( "CD4 Tcell" ).type;

        // (Adrianne) if there is a T cell in a mac's neighbourhood AND a mac has already begin phagocytosing, then there will be changes to the macs actions 
        //        int n = 0; 
        //        Cell pContactCell = neighbors[n]; 
        //        while( n < neighbors.size() )
        //        {
        //            pContactCell = neighbors[n]; 
        for( Cell pContactCell : neighbors )
        {
            double cell_cell_distance = Math
                    .sqrt( ( pContactCell.position[0] - pCell.position[0] ) * ( pContactCell.position[0] - pCell.position[0] )
                            + ( pContactCell.position[1] - pCell.position[1] ) * ( pContactCell.position[1] - pCell.position[1] ) );
            double radius_mac = pCell.phenotype.geometry.radius; // (Adrianne) radius of DC)
            double radius_test_cell = pContactCell.phenotype.geometry.radius; // (Adrianne) radius of test cell)

            // (Adrianne) if it is not me, not dead and is a CD8 T cell that is within a very short distance from me, I will stop secreting pro-inflammatory cytokine
            if( pContactCell != pCell && pContactCell.phenotype.death.dead == false && pContactCell.type == CD8_Tcell_type
                    && pCell.customData.get( "activated_immune_cell" ) > 0.5
                    && cell_cell_distance <= pCell.getModel().getParameterDouble( "epsilon_distance" )
                            * ( radius_mac + radius_test_cell ) )
            {
                phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = 0;// Contact with CD8 T cell turns off pro-inflammatory cytokine secretion
                phenotype.secretion.secretionRates[antiinflammatory_cytokine_index] = pCell.customData
                        .get( "antiinflammatory_cytokine_secretion_rate_by_macrophage" );// and turns on anti-inflammatory cytokine secretion
                //                n=neighbors.size();
                break;
            }
            // (Adrianne) if it is not me, not dead and is a CD4 T cell that is within a very short distance from me, I will be able to phagocytose infected (but not neccesarily dead) cells
            else if( pContactCell != pCell && pContactCell.phenotype.death.dead == false && pContactCell.type == CD4_Tcell_type
                    && pCell.customData.get( "activated_immune_cell" ) > 0.5
                    && cell_cell_distance <= pCell.getModel().getParameterDouble( "epsilon_distance" )
                            * ( radius_mac + radius_test_cell ) )
            {
                pCell.customData.set( "ability_to_phagocytose_infected_cell", 1 ); // (Adrianne) contact with CD4 T cell induces macrophage's ability to phagocytose infected cells
                //                n=neighbors.size();
                break;
            }
            //            n++;
        }

        // (Adrianne) if macrophage volume exceeds a threshold value we say it is "exhausted" and unable to phagocytose until it's volume drops below this threshold
        if( pCell.phenotype.volume.total > pCell.customData.get( "threshold_macrophage_volume" ) )
        {
            // (Adrianne) when a macrophage is in an exhausted state it has a death rate  2.1e-4
            phenotype.death.rates.set( apoptosis_index, pCell.customData.get( "exhausted_macrophage_death_rate" ) );
            return;
        }

        // (Adrianne) obtain index for tracking time when next phagocytosis event is possible
        int time_to_next_phagocytosis_index = pCell.customData.findVariableIndex( "time_to_next_phagocytosis" );
        // (Adrianne) check if still phagocytosing something, added if statement to say that if cell is still internalising current material not to phagocytose anything else
        if( pCell.customData.variables.get( time_to_next_phagocytosis_index ).value > pCell.getModel().getCurrentTime() )
        {
            return;
        }

        double probability_of_phagocytosis = pCell.customData.get( "phagocytosis_rate" ) * dt;
        /* // remove in v 3.2 
        double max_phagocytosis_volume = pCell.customData.get("phagocytosis_relative_target_cutoff_size" ] * pCD.phenotype.volume.total; 
         */
        // (Adrianne) add an additional variable that is the time taken to ingest material 
        double material_internalisation_rate = pCell.customData.get( "material_internalisation_rate" );

        //            n = 0; 
        //            Cell pTestCell = neighbors[n]; 
        //            while( n < neighbors.size() )
        //            {
        //                pTestCell = neighbors[n]; 
        for( Cell pTestCell : neighbors )
        {
            int nP = pTestCell.customData.findVariableIndex( "viral_protein" ); //(Adrianne) finding the viral protein inside cells
            // if it is not me and not a macrophage 
            if( pTestCell != pCell && pTestCell.phenotype.death.dead == true
                    && pCell.getModel().getRNG().UniformRandom() < probability_of_phagocytosis ) // && // remove in v 3.2 
            //          pTestCell.phenotype.volume.total < max_phagocytosis_volume ) / remove in v 3.2 
            {
                {
                    // (Adrianne) obtain volume of cell to be ingested
                    double volume_ingested_cell = pTestCell.phenotype.volume.total;

                    pCell.ingestCell( pTestCell );

                    // (Adrianne)(assume neutrophils same as macrophages) neutrophils phagocytose material 1micron3/s so macrophage cannot phagocytose again until it has elapsed the time taken to phagocytose the material
                    double time_to_ingest = volume_ingested_cell * material_internalisation_rate;// convert volume to time taken to phagocytose
                    // (Adrianne) update internal time vector in macrophages that tracks time it will spend phagocytosing the material so they can't phagocytose again until this time has elapsed
                    pCell.customData.variables.get( time_to_next_phagocytosis_index ).value = pCell.getModel().getCurrentTime()
                            + time_to_ingest;
                }

                // activate the cell 
                phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = pCell.customData
                        .get( "activated_cytokine_secretion_rate" ); // 10;
                phenotype.secretion.saturationDensities[proinflammatory_cytokine_index] = 1;

                phenotype.secretion.uptakeRates[proinflammatory_cytokine_index] = 0.0;

                phenotype.motility.migrationSpeed = pCell.customData.get( "activated_speed" );

                //(adrianne V5) adding virus uptake by phagocytes
                phenotype.secretion.uptakeRates[virus_index] = pCell.getModel().getParameterDouble( "phagocytes_virus_uptake_rate" );

                pCell.customData.set( "activated_immune_cell", 1.0 );
//                System.out.println( "Activated dead" );

                return;
            }
            else if( pTestCell != pCell && pCell.customData.get( "ability_to_phagocytose_infected_cell" ) == 1
                    && pTestCell.customData.get( nP ) > 1 && pCell.getModel().getRNG().UniformRandom() < probability_of_phagocytosis ) // (Adrianne) macrophages that have been activated by T cells can phagocytose infected cells that contain at least 1 viral protein
            {
                {
                    // (Adrianne) obtain volume of cell to be ingested
                    double volume_ingested_cell = pTestCell.phenotype.volume.total;

                    pCell.ingestCell( pTestCell );

                    // (Adrianne)(assume neutrophils same as macrophages) neutrophils phagocytose material 1micron3/s so macrophage cannot phagocytose again until it has elapsed the time taken to phagocytose the material
                    double time_to_ingest = volume_ingested_cell * material_internalisation_rate;// convert volume to time taken to phagocytose
                    // (Adrianne) update internal time vector in macrophages that tracks time it will spend phagocytosing the material so they can't phagocytose again until this time has elapsed
                    pCell.customData.variables.get( time_to_next_phagocytosis_index ).value = pCell.getModel().getCurrentTime()
                            + time_to_ingest;
                }

                // activate the cell 
                phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = pCell.customData
                        .get( "activated_cytokine_secretion_rate" ); // 10;
                phenotype.secretion.saturationDensities[proinflammatory_cytokine_index] = 1;

                phenotype.secretion.uptakeRates[proinflammatory_cytokine_index] = 0.0;

                //(adrianne v5) adding virus uptake by phagocytes
                phenotype.secretion.uptakeRates[virus_index] = pCell.getModel().getParameterDouble( "phagocytes_virus_uptake_rate" );

                phenotype.motility.migrationSpeed = pCell.customData.get( "activated_speed" );

                pCell.customData.set( "activated_immune_cell", 1.0 );
//                System.out.println( "Activated alice" );

                return;
            }
            else if( pTestCell != pCell && pTestCell.phenotype.death.dead == false
                    && pTestCell.phenotype.molecular.internSubstrates[antibody_index] > 1e-12 ) //  .internalized_total_substrates[antibody_index]>1e-12 ) 
            // (Adrianne V5) macrophages can phaogyctose infected cell if it has some non-trivial bound antibody
            {
                double antibody_level = pTestCell.phenotype.molecular.internSubstrates[antibody_index];
                double antibody_half_effect = pCell.getModel().getParameterDouble( "antibody_half_effect" );
                // (AJ-V5) recalculate probability of phagocytosis based on bound anitbody
                if( pCell.getModel().getRNG().UniformRandom() < probability_of_phagocytosis
                        + ( 1 - probability_of_phagocytosis ) * ( antibody_level ) / ( antibody_level + antibody_half_effect ) )
                // cell phagocytses
                {
                    // (Adrianne) obtain volume of cell to be ingested
                    double volume_ingested_cell = pTestCell.phenotype.volume.total;

                    pCell.ingestCell( pTestCell );

                    // (Adrianne)(assume neutrophils same as macrophages) neutrophils phagocytose material 1micron3/s so macrophage cannot phagocytose again until it has elapsed the time taken to phagocytose the material
                    double time_to_ingest = volume_ingested_cell * material_internalisation_rate;// convert volume to time taken to phagocytose
                    // (Adrianne) update internal time vector in macrophages that tracks time it will spend phagocytosing the material so they can't phagocytose again until this time has elapsed
                    pCell.customData.variables.get( time_to_next_phagocytosis_index ).value = pCell.getModel().getCurrentTime()
                            + time_to_ingest;
                }
            }
            //                n++; 
        }
        return;
    }
}