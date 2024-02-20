package ru.biosoft.physicell.sample_projects.mechano;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class CancerPhenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        if( SignalBehavior.getSingleSignal( pCell, "dead" ) > 0.5 )
        {
            return;
        }

        double b = SignalBehavior.get_single_base_behavior( pCell, "cycle entry" );
        if( SignalBehavior.getSingleSignal( pCell, "pressure" ) > 0.75 )
        {
            b = 0;
        }
        SignalBehavior.setSingleBehavior( pCell, "cycle entry", b );

        if( SignalBehavior.getSingleSignal( pCell, "time" ) > 10000 )
        {
            SignalBehavior.setSingleBehavior( pCell, "apoptosis", 9e99 );
            SignalBehavior.setSingleBehavior( pCell, "cell detachment rate", 9e9 );
        }
    }
}