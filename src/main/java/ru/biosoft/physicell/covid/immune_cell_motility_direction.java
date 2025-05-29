package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;

public class immune_cell_motility_direction extends UpdateMigrationBias
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        //    void immune_cell_motility_direction(Cell pCell, Phenotype phenotype, double dt)
        //    {
        if( phenotype.death.dead == true )
        {
            phenotype.motility.migrationSpeed = 0.0;
            return;
        }

        int chemokine_index = pCell.getMicroenvironment().findDensityIndex( "chemokine" );
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        // if not activated, chemotaxis along debris 

        phenotype.motility.migrationBiasDirection = pCell.nearest_gradient( debris_index );
        VectorUtil.normalize( phenotype.motility.migrationBiasDirection );
        if( pCell.customData.get( "activated_immune_cell" ) < 0.5 )
        {
            return;
        }

        // if activated, follow the weighted direction 

        VectorUtil.prod( phenotype.motility.migrationBiasDirection, pCell.customData.get( "sensitivity_to_debris_chemotaxis" ) );
        //        phenotype.motility.migrationBiasDirection *= pCell.customData.get("sensitivity_to_debris_chemotaxis");

        double[] gradC = pCell.nearest_gradient( chemokine_index );
        VectorUtil.normalize( gradC );
        VectorUtil.prod( gradC, pCell.customData.get( "sensitivity_to_chemokine_chemotaxis" ) );
        //        gradC *= pCell.customData.get("sensitivity_to_chemokine_chemotaxis"];

        VectorUtil.sum( phenotype.motility.migrationBiasDirection, gradC );
        //        phenotype.motility.migrationBiasDirection += gradC; 

        VectorUtil.normalize( phenotype.motility.migrationBiasDirection );
    }
}