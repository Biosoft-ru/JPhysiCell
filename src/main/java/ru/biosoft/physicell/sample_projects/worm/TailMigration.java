package ru.biosoft.physicell.sample_projects.worm;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.standard.Chemotaxis;

public class TailMigration extends Chemotaxis
{
    private int direction;
    private double speed;
    private double bias;
    private double persistenceTime;

    public TailMigration(Model model)
    {
        direction = model.getParameterInt( "tail_migration_direction" );
        speed = model.getParameterDouble( "tail_migration_speed" );
        bias = model.getParameterDouble( "tail_migration_bias" );
        persistenceTime = model.getParameterDouble( "tail_migration_persistence" );
    }

    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        phenotype.motility.chemotaxisDirection = direction;
        phenotype.motility.migrationSpeed = speed;
        phenotype.motility.migrationBias = bias;
        phenotype.motility.persistenceTime = persistenceTime;
        super.execute( pCell, phenotype, dt );
    }
}