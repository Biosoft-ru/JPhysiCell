package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;

public class WorkerCellRule extends CustomCellRule
{
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        // if I am dead, don't bother
        if( SignalBehavior.getSingleSignal( pCell, "dead" ) > 0.5 )
        {
            // the cell death functions don't automatically turn off custom functions,
            // since those are part of mechanics.

            // Let's just fully disable now.
            pCell.functions.customCellRule = null;
            return;
        }

        // am I searching for cargo? if so, see if I've found it
        if( pCell.state.numberAttachedCells() == 0 )
        {
            boolean attached = false; // want to limit to one attachment
            for( Cell nearbyCell : pCell.cells_in_my_container() )
            {
                if( nearbyCell == pCell )
                    continue;
                // if it is expressing the receptor, dock with it
                if( SignalBehavior.getSingleSignal( nearbyCell, "custom:receptor" ) > 0.5 && attached == false )
                {
                    Cell.attachcCells( pCell, nearbyCell );
                    // nearby[i].custom_data["receptor"] = 0.0; // put into cargo cell rule instead?
                    // nearby[i].phenotype.secretion.set_all_secretion_to_zero(); // put into cargo rule instead?
                    attached = true;
                    break;
                }
            }
        }

        // from prior motility function
        double o2 = SignalBehavior.getSingleSignal( pCell, "oxygen" );
        double chemoattractant = SignalBehavior.getSingleSignal( pCell, "chemoattractant" );
        double detectionThreshold = SignalBehavior.getSingleSignal( pCell, "custom:motility_shutdown_detection_threshold" );

        // if attached, biased motility towards director chemoattractant
        // otherwise, biased motility towards cargo chemoattractant
        double attachedMigrationBias = SignalBehavior.getSingleSignal( pCell, "custom:attached_worker_migration_bias" );
        double unattachedMigrationBias = SignalBehavior.getSingleSignal( pCell, "custom:unattached_worker_migration_bias" );

        if( pCell.state.numberAttachedCells() > 0 )
        {
            SignalBehavior.setSingleBehavior( pCell, "migration bias", attachedMigrationBias );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to oxygen", -1 );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to chemoattractant", 0 );
        }
        else
        {
            // if there is no detectable signal, shut down motility (permanently)
            if( chemoattractant < detectionThreshold )
            {
                SignalBehavior.setSingleBehavior( pCell, "migration speed", 0 );
            }
            SignalBehavior.setSingleBehavior( pCell, "migration bias", unattachedMigrationBias );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to oxygen", 0 );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to chemoattractant", 1 );
        }
    }
}