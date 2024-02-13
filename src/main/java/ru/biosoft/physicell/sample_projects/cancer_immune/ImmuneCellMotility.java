package ru.biosoft.physicell.sample_projects.cancer_immune;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.update_migration_bias;

public class ImmuneCellMotility implements update_migration_bias
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        // if attached, biased motility towards director chemoattractant 
        // otherwise, biased motility towards cargo chemoattractant 
        int immuneFactorIndex = microenvironment.findDensityIndex( "immunostimulatory factor" );

        // if not docked, attempt biased chemotaxis 
        if( pCell.state.attachedCells.size() == 0 )
        {
            phenotype.motility.isMotile = true;
            phenotype.motility.migrationBiasDirection = pCell.nearest_gradient( immuneFactorIndex ).clone();
            VectorUtil.normalize( ( phenotype.motility.migrationBiasDirection ) );
        }
        else
        {
            phenotype.motility.isMotile = false;
        }
    }
}