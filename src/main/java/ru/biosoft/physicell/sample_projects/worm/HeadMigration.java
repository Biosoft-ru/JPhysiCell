package ru.biosoft.physicell.sample_projects.worm;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.standard.Chemotaxis;

public class HeadMigration extends Chemotaxis
{
    private int direction;
    private double speed;
    private double bias;
    private double persistenceTime;

    public HeadMigration(Model model)
    {
        direction = model.getParameterInt( "head_migration_direction" );
        speed = model.getParameterDouble( "head_migration_speed" );
        bias = model.getParameterDouble( "head_migration_bias" );
        persistenceTime = model.getParameterDouble( "head_migration_persistence" );
    }

    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        phenotype.motility.chemotaxisDirection = direction;
        phenotype.motility.migrationSpeed = speed;
        phenotype.motility.migrationBias = bias;
        phenotype.motility.persistenceTime = persistenceTime;
        // use this for fun rotational paths 
        /*
        double r = norm( pCell.position ) + 1e-16; 
        phenotype.motility.migration_bias_direction[0] = - pCell.position[1] / r; 
        phenotype.motility.migration_bias_direction[1] = pCell.position[0] / r; 
        
        normalize( &(phenotype.motility.migration_bias_direction) ); 
        return; 
        */
        super.execute( pCell, phenotype, dt );
    }
}