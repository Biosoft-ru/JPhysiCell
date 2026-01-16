package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;


public class CD4_Tcell_mechanics extends CustomCellRule
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
        if( Util.check_for_out_of_bounds( pCell, 10.0 ) )
        {
            //		#pragma omp critical
            ( (MicrometsModel)pCell.getModel() ).addCellsToMoveFromEdge( pCell );
        }

        return;
    }
}
