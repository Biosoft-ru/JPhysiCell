package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;

public class AdvancedChemotaxisNormalized implements UpdateMigrationBias
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // We'll work directly on the migration bias direction 
        double[] pVec = phenotype.motility.migrationBiasDirection;
        // reset to zero. use memset to be faster??
        for( int i = 0; i < 3; i++ )
            pVec[i] = 0;

        // a place to put each gradient prior to normalizing it 
        double[] temp = new double[3];
        // weighted combination of the gradients 
        for( int i = 0; i < phenotype.motility.chemotacticSensitivities.length; i++ )
        {
            // get and normalize ith gradient 
            temp = pCell.nearest_gradient( i );
            VectorUtil.normalize( temp );
            VectorUtil.axpy( pVec, phenotype.motility.chemotacticSensitivities[i], temp );
        }
        // normalize that 
        VectorUtil.normalize( pVec );
    }
}