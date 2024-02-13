package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class CargoPhenotype extends UpdatePhenotype
{
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        // if dettached and receptor on, secrete signal
        // if dettached and receptor off, secrete chemo
        double receptor = SignalBehavior.getSingleSignal( pCell, "custom:receptor" );

        if( pCell.state.numberAttachedCells() == 0 )
        {
            if( receptor > 0.1 )
            {
                SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 10 );
                SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 0 );
            }
            else
            {
                SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 0 );
                SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 10 );
            }
            return;
        }

        // if you reach this point of the code, the cell is attached
        // if attached and oxygen high, secrete nothing, receptor off
        // if attached and oxygen low, dettach, start secreting chemo, receptor off
        double o2 = SignalBehavior.getSingleSignal( pCell, "oxygen" );
        double o2_drop = SignalBehavior.getSingleSignal( pCell, "custom:cargo_release_o2_threshold" );

        if( o2 > o2_drop )
        {
            SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 0 );
            SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 0 );
            SignalBehavior.setSingleBehavior( pCell, "custom:receptor", 0 );
        }
        else
        {
            SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 0 );
            SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 10 );
            SignalBehavior.setSingleBehavior( pCell, "custom:receptor", 0 );
            pCell.removeAllAttachedCells();
        }
    }
}