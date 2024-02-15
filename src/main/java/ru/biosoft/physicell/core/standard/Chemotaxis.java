package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;

public class Chemotaxis implements UpdateMigrationBias
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // bias direction is gradient for the indicated substrate 
        phenotype.motility.migrationBiasDirection = pCell.nearest_gradient( phenotype.motility.chemotaxisIndex ).clone();
        // move up or down gradient based on this direction 
        VectorUtil.prod( phenotype.motility.migrationBiasDirection, phenotype.motility.chemotaxisDirection );
        VectorUtil.normalize( phenotype.motility.migrationBiasDirection );
    }
}