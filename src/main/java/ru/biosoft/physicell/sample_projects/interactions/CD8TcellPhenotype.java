package ru.biosoft.physicell.sample_projects.interactions;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.BasicSignaling;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;

public class CD8TcellPhenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( pCell.typeName );
        Microenvironment microenvironment = pCell.getMicroenvironment();
        //        int nR = microenvironment.findDensityIndex( "resource" );
        //        int nTox = microenvironment.findDensityIndex( "toxin" );
        int nDebris = microenvironment.findDensityIndex( "debris" );
        int nPIF = microenvironment.findDensityIndex( "pro-inflammatory" );

        double[] samples = pCell.nearest_density_vector();
        double PIF = samples[nPIF];

        // if dead, release debris
        if( phenotype.death.dead == true )
        {
            phenotype.secretion.netExportRates[nDebris] = phenotype.volume.total;
            pCell.functions.updatePhenotype = null;
            return;
        }

        // migration bias increases with pro-inflammatory 
        //        double signal = PIF;
        double base_val = pCD.phenotype.motility.migrationBias;
        double max_val = 0.75;
        double half_max = pCD.custom_data.get( "migration_bias_halfmax" ); // 0.05 // 0.25 
        double hill = BasicSignaling.Hill_response_function( PIF, half_max, 1.5 );
        phenotype.motility.migrationBias = base_val + ( max_val - base_val ) * hill;
    }

    @Override
    public String display()
    {
        return "Pro-inflammatory increases migration." + " Release debris upon death.";
    }
}