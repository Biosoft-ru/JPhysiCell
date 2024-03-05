package ru.biosoft.physicell.sample_projects.worm;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;

public class MiddleMigration extends UpdateMigrationBias
{
    double speed;

    public MiddleMigration(Model model)
    {
        speed = model.getParameterDouble( "middle_migration_speed" );
    }

    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // get velocity from "Upstream" 
        Cell headCell = null;
        for( Cell cell : pCell.state.attachedCells )
        {
            if( headCell == null || cell.customData.get( "head" ) > headCell.customData.get( "head" ) )
                headCell = cell;
        }
        phenotype.motility.migrationSpeed = speed;
        phenotype.motility.migrationBiasDirection = headCell.phenotype.motility.migrationBiasDirection;
        VectorUtil.normalize( phenotype.motility.migrationBiasDirection );
    }
}