package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import java.util.Set;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.CellFunctions.custom_cell_rule;

public class WorkerCellRule implements custom_cell_rule
{
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        // if I am dead, don't bother
        if( SignalBehavior.getSingleSignal( pCell, "dead" ) > 0.5 )
        {
            // the cell death functions don't automatically turn off custom functions,
            // since those are part of mechanics.

            // Let's just fully disable now.
            pCell.functions.custom_cell_rule = null;
            return;
        }

        // am I searching for cargo? if so, see if I've found it
        if( pCell.state.numberAttachedCells() == 0 )
        {
            Set<Cell> nearby = pCell.cells_in_my_container();
            boolean attached = false; // want to limit to one attachment
            for( Cell nearbyCell : nearby )
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
        double detection_threshold = SignalBehavior.getSingleSignal( pCell, "custom:motility_shutdown_detection_threshold" );

        // if attached, biased motility towards director chemoattractant
        // otherwise, biased motility towards cargo chemoattractant
        double attached_worker_migration_bias = SignalBehavior.getSingleSignal( pCell, "custom:attached_worker_migration_bias" );
        double unattached_worker_migration_bias = SignalBehavior.getSingleSignal( pCell, "custom:unattached_worker_migration_bias" );

        if( pCell.state.numberAttachedCells() > 0 )
        {
            SignalBehavior.setSingleBehavior( pCell, "migration bias", attached_worker_migration_bias );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to oxygen", -1 );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to chemoattractant", 0 );
        }
        else
        {
            // if there is no detectable signal, shut down motility (permanently)
            if( chemoattractant < detection_threshold )
            {
                SignalBehavior.setSingleBehavior( pCell, "migration speed", 0 );
            }
            SignalBehavior.setSingleBehavior( pCell, "migration bias", unattached_worker_migration_bias );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to oxygen", 0 );
            SignalBehavior.setSingleBehavior( pCell, "chemotactic response to chemoattractant", 1 );
        }
    }
}