package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class CD8_Tcell_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        // MISSING CHANGE DEATH RATE (DEPEND ON GENERATION)
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris");

        if( phenotype.death.dead == true )
        {
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get("debris_secretion_rate");
            return;
        }
        return;
    }
}
