package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;
import ru.biosoft.physicell.core.Phenotype;

public class AdvancedChemotaxis extends UpdateMigrationBias
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        double[] pVec = phenotype.motility.migrationBiasDirection;
        VectorUtil.zero( pVec );
        for( int i = 0; i < phenotype.motility.chemotacticSensitivities.length; i++ )
            VectorUtil.axpy( pVec, phenotype.motility.chemotacticSensitivities[i], pCell.nearest_gradient( i ) );
        VectorUtil.normalize( pVec );
    }

    @Override
    public String display()
    {
        return "Advanced chemotaxis (weighted combintaion of gradients)";
    }
}