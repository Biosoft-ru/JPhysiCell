package ru.biosoft.physicell.sample_projects.heterogeneity;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.StandardModels.O2based;

public class TumorPhenotype extends O2based
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        super.execute( pCell, phenotype, dt );

        if( SignalBehavior.getSingleSignal( pCell, "dead" ) > 0.5 )
        {
            pCell.functions.updatePhenotype = null;// if cell is dead, don't bother with future phenotype changes. 
        }
        else
        {
            // multiply proliferation rate by the oncoprotein 
            double rate = SignalBehavior.getSinglBehavior( pCell, "cycle entry" );
            double factor = SignalBehavior.getSingleSignal( pCell, "custom:oncoprotein" );
            SignalBehavior.setSingleBehavior( pCell, "cycle entry", rate * factor );
        }
    }
}