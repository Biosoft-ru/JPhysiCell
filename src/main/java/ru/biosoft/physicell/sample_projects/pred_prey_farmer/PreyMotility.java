package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;
import ru.biosoft.physicell.core.Phenotype;

public class PreyMotility implements UpdateMigrationBias
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        weighted_motility_function( pCell, phenotype, dt );
    }

    public static void weighted_motility_function(Cell pCell, Phenotype phenotype, double dt)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        // find the indices for each major substrate 
        int prey_index = microenvironment.findDensityIndex( "prey signal" );
        int predator_index = microenvironment.findDensityIndex( "predator signal" );
        int food_index = microenvironment.findDensityIndex( "food" );

        // zero out the motility bias direction. use a pointer to make this easier 
        double[] pV = phenotype.motility.migrationBiasDirection;
        pV[0] = 0;
        pV[1] = 0;
        pV[2] = 0;
        //  (pV) = {0,0,0}; // pCell.position; 
        //*pV *= -0.00001; 

        // v += prey_weight * grad(prey) 
        VectorUtil.axpy( pV, pCell.custom_data.get( "prey_weight" ), pCell.nearest_gradient( prey_index ) );
        // v += predator_weight * grad(predator) 
        VectorUtil.axpy( pV, pCell.custom_data.get( "predator_weight" ), pCell.nearest_gradient( predator_index ) );
        // v += food_weight * grad(food) 
        VectorUtil.axpy( pV, pCell.custom_data.get( "food_weight" ), pCell.nearest_gradient( food_index ) );

        VectorUtil.normalize( pV );
    }
}