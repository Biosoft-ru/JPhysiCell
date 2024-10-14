package ru.biosoft.physicell.sample_projects.mechano;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;

public class CancerPhenotype extends UpdatePhenotype
{
    SignalBehavior signals;

    public CancerPhenotype(Model model)
    {
        signals = model.getSignals();
    }

    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        if( signals.getSingleSignal( pCell, "dead" ) > 0.5 )
        {
            return;
        }

        double b = signals.getSinglBaseBehavior( pCell, "cycle entry" );
        if( signals.getSingleSignal( pCell, "pressure" ) > 0.75 )
        {
            b = 0;
        }
        signals.setSingleBehavior( pCell, "cycle entry", b );

        if( signals.getSingleSignal( pCell, "time" ) > 10000 )
        {
            signals.setSingleBehavior( pCell, "apoptosis", 9e99 );
            signals.setSingleBehavior( pCell, "cell detachment rate", 9e9 );
        }
    }
}