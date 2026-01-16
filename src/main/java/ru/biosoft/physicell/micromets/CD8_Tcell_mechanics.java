package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.RandomGenerator;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;

public class CD8_Tcell_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {

        RandomGenerator rng = pCell.getModel().getRNG();
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

        // if I am not adhered to a cell, turn motility on
        if( pCell.state.neighbors.size() == 0 )
        {
            phenotype.motility.isMotile = true;
        }
        else
        {
            phenotype.motility.isMotile = false;
        }

        // if I'm adhered to something ...
        if( pCell.state.numberAttachedCells() > 0 )
        {
            // decide whether to detach
            boolean detach_me = false;

            if( rng.UniformRandom() < dt / ( pCell.customData.get( "cell_attachment_lifetime" ) + 1e-15 ) )
            {
                detach_me = true;
            }

            // if I detach, go through the process
            if( detach_me )
            {
                pCell.removeAllAttachedCells();
                // resume motile behavior
                phenotype.motility.isMotile = true;
            }
            return;
        }

        // I'm not attached, look for cells nearby and try to attach
        // if this returns non-NULL, we're now attached to a cell
        if( Util.immune_cell_check_neighbors_for_attachment( pCell, dt ) != null )
        {
            // set motility off
            phenotype.motility.isMotile = false;
            return;
        }

        return;
    }

    


}