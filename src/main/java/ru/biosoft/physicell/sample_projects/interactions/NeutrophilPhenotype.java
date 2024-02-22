package ru.biosoft.physicell.sample_projects.interactions;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.BasicSignaling;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;

/**
 * Migrate to pro-inflammatory, when die - release debris
 */
public class NeutrophilPhenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( pCell.typeName );
        Microenvironment m = pCell.getMicroenvironment();

        if( phenotype.death.dead == true )
        {
            int nDebris = m.findDensityIndex( "debris" );
            phenotype.secretion.netExportRates[nDebris] = phenotype.volume.total;
            pCell.functions.updatePhenotype = null;
            return;
        }

        // migration bias increases with pro-inflammatory
        int nPIF = m.findDensityIndex( "pro-inflammatory" );
        double PIF = pCell.nearest_density_vector()[nPIF];
        double base_val = pCD.phenotype.motility.migrationBias;
        double max_val = 0.75;
        double half_max = pCD.custom_data.get( "migration_bias_halfmax" ); // 0.25 
        double hill = BasicSignaling.Hill_response_function( PIF, half_max, 1.5 );
        phenotype.motility.migrationBias = base_val + ( max_val - base_val ) * hill;
    }

    @Override
    public String display()
    {
        return "Pro-inflammatory increases migration bias. " + "Debris are released upon death.";
    }
}