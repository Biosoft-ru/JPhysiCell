package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.Phenotype;

public class DC_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
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
            ( (MicrometsModel)pCell.getModel() ).addCellsToMoveFromEdge( pCell );
        }

        return;
    }
}
