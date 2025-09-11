package ru.biosoft.physicell.covid;

import java.util.Set;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class neutrophil_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        int apoptosis_index = phenotype.death.findDeathModelIndex( "apoptosis" );
        CellDefinition pCD = pCell.getModel().getCellDefinition( "neutrophil" );
        int proinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "pro-inflammatory cytokine" );
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );
        int chemokine_index = pCell.getMicroenvironment().findDensityIndex( "chemokine" );
        int ROS_index = pCell.getMicroenvironment().findDensityIndex( "ROS" );
        int virus_index = pCell.getMicroenvironment().findDensityIndex( "virion" );

        if( phenotype.death.dead == true )
        {
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
            return;
        }

        // check for cells to eat 
        Set<Cell> neighbors = pCell.cells_in_my_container();

        // at least one of the cells is pCell 
        if( neighbors.size() < 2 )
            return;

        // (Adrianne) if neutrophil volume exceeds a threshold value we say it is "exhausted" and unable to phagocytose until it's volume drops below this threshold
        if( pCell.phenotype.volume.total > pCell.customData.get( "threshold_neutrophil_volume" ) )
            return;

        // (Adrianne) obtain index for tracking time to next phagocytosis event is possible
        int time_to_next_phagocytosis_index = pCell.customData.findVariableIndex( "time_to_next_phagocytosis" );
        // (Adrianne) check if still phagocytosing something, added if statement to say that if cell is still internalising current material not to phagocytose anything else
        if( pCell.customData.variables.get( time_to_next_phagocytosis_index ).value > pCell.getModel().getCurrentTime() )
            return;

        double probability_of_phagocytosis = pCell.customData.get( "phagocytosis_rate" ) * dt;
        double max_phagocytosis_volume = pCell.customData.get( "phagocytosis_relative_target_cutoff_size" ) * pCD.phenotype.volume.total;

        // (Adrianne) add an additional variable that is the time taken to ingest material 
        double material_internalisation_rate = pCell.customData.get( "material_internalisation_rate" );

        for( Cell pTestCell : neighbors )
        {
            // if it is not me and the target is dead 
            if( pTestCell != pCell && pTestCell.phenotype.death.dead
                    && pCell.getModel().getRNG().UniformRandom() < probability_of_phagocytosis
                    && pTestCell.phenotype.volume.total < max_phagocytosis_volume )
            {

                // (Adrianne) obtain volume of cell to be ingested
                double volume_ingested_cell = pTestCell.phenotype.volume.total;

                // remove_all_adhesions( pTestCell ); // debug 
                pCell.ingestCell( pTestCell );

                // (Adrianne)(assume neutrophils same as macrophages) neutrophils phagocytose material 1micron3/s so macrophage cannot phagocytose again until it has elapsed the time taken to phagocytose the material
                double time_to_ingest = volume_ingested_cell * material_internalisation_rate;// convert volume to time taken to phagocytose
                // (Adrianne) update internal time vector in macrophages that tracks time it will spend phagocytosing the material so they can't phagocytose again until this time has elapsed
                pCell.customData.variables.get( time_to_next_phagocytosis_index ).value = pCell.getModel().getCurrentTime()
                        + time_to_ingest;

                // activate the cell 
                phenotype.secretion.secretionRates[proinflammatory_cytokine_index] = pCell.customData
                        .get( "activated_cytokine_secretion_rate" ); // 10;
                phenotype.secretion.saturationDensities[proinflammatory_cytokine_index] = 1;

                // (Adrianne V5) Cell starts secreting ROS
                phenotype.secretion.secretionRates[ROS_index] = pCell.getModel().getParameterDouble( "ROS_secretion_rate" ); // 10;
                phenotype.secretion.saturationDensities[proinflammatory_cytokine_index] = 1;

                //(adrianne V5) adding virus uptake by phagocytes
                phenotype.secretion.uptakeRates[virus_index] = pCell.getModel().getParameterDouble( "phagocytes_virus_uptake_rate" );

                phenotype.motility.migrationSpeed = pCell.customData.get( "activated_speed" );

                pCell.customData.set( "activated_immune_cell", 1.0 );
                return;
            }
        }
    }
}