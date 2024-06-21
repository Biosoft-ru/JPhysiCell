package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Model.Event;
import ru.biosoft.physicell.core.PhysiCellUtilities;

public class TherapyEvent extends Event
{
    public TherapyEvent(double executionTime)
    {
        super( executionTime );
    }

    @Override
    public void execute(Model model) throws Exception
    {
        //        System.out.println( "Therapy started!" );
        model.setSaveFullInterval( model.getParameterDouble( "save_interval_after_therapy_start" ) ); // 3.0; 
        introduceBiorobots( model );
    }

    public void introduceBiorobots(Model model) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        // idea: we'll "inject" them in a little column
        double workerFraction = model.getParameterDouble( "worker_fraction" ); // 0.10; /* param */
        int numberInjectedCells = model.getParameterInt( "number_of_injected_cells" ); // 500; /* param */

        // make these vary with domain size
        double left = m.options.X_range[1] - 150.0; // 600.0;
        double right = m.options.X_range[1] - 50.0; // 700.0;
        double bottom = m.options.Y_range[0] + 50.0; // -700;
        double top = m.options.Y_range[1] - 50.0; // 700;

        CellDefinition workerCD = CellDefinition.getCellDefinition( "worker cell" );
        CellDefinition cargoCD = CellDefinition.getCellDefinition( "cargo cell" );

        for( int i = 0; i < numberInjectedCells; i++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( left, right );
            position[1] = PhysiCellUtilities.UniformRandom( bottom, top );

            if( PhysiCellUtilities.UniformRandom() <= workerFraction )
            {
                Cell.createCell( workerCD, m, position );
            }
            else
            {
                Cell.createCell( cargoCD, m, position );
            }
        }
    }
}