package ru.biosoft.physicell.sample_projects.heterogeneity;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.StandardModels.update_cell_and_death_parameters_O2_based;

public class TumorPhenotype extends update_cell_and_death_parameters_O2_based
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        super.execute( pCell, phenotype, dt );

        if( SignalBehavior.get_single_signal( pCell, "dead" ) > 0.5 )
        {
            pCell.functions.updatePhenotype = null;// if cell is dead, don't bother with future phenotype changes. 
        }
        else
        {
            // multiply proliferation rate by the oncoprotein 
            double rate = SignalBehavior.get_single_behavior( pCell, "cycle entry" );
            double factor = SignalBehavior.get_single_signal( pCell, "custom:oncoprotein" );
            //            System.out.println( pCell.toString() + " " + rate + " " + factor + " " + rate * factor );
            //            System.out
            //                    .println( pCell.toString() + " " + cycle_rate + " " + SignalBehavior.get_single_signal( pCell, "custom:oncoprotein" ) );
            rate *= SignalBehavior.get_single_signal( pCell, "custom:oncoprotein" );
            //            System.out.println( pCell.toString() + "\t" + rate );
            SignalBehavior.setSingleBehavior( pCell, "cycle entry", rate );
        }
    }

    //    private void write(String path)
    //    {
    //        File f= new File(path);
    //        try(BufferedWriter bw = new B)
    //    }
}