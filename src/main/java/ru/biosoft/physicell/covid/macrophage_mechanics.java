package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;

public class macrophage_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        if( phenotype.death.dead == true )
        {
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
            return;
        }

        // bounds check 
        if( ModelCovid.check_for_out_of_bounds( pCell, 10.0 ) )
        {
            //            #pragma omp critical 
            {
                ( (ModelCovid)pCell.getModel() ).cells_to_move_from_edge.add( pCell );
            }
            // replace_out_of_bounds_cell( pCell, 10.0 );
            // return; 
        }

        //  // death check 
        //  if( phenotype.death.dead == true ) 
        //  { remove_all_adhesions( pCell ); }

        //        return; 
    }
}