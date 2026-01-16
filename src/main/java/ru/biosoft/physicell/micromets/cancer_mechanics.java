package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;


public class cancer_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        // if I'm dead, don't bother
        if( phenotype.death.dead == true )
        {
            // the cell death functions don't automatically turn off custom functions,
            // since those are part of mechanics.
            // remove_all_adhesions( pCell );

            // Let's just fully disable now.
            pCell.functions.customCellRule = null;
            pCell.functions.contact = null;

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
            return;
        }

        // bounds check
        if( Util.check_for_out_of_bounds( pCell, 10.0 ) )
        {
            ( (MicrometsModel)pCell.getModel() ).addCellsToMoveFromEdge( pCell );
            //            #pragma omp critical
            //            {
            //                cells_to_move_from_edge.push_back( pCell );
            //            }
        }

        return;
    }
}
