package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class fibroblast_phenotype extends UpdatePhenotype
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        int antiinflammatory_cytokine_index = pCell.getMicroenvironment().findDensityIndex( "anti-inflammatory cytokine" );
        int collagen_index = pCell.getMicroenvironment().findDensityIndex( "collagen" );
        double TGF_beta = pCell.nearest_density_vector()[antiinflammatory_cytokine_index];
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );
        CellDefinition pCD = pCell.getModel().getCellDefinition( "fibroblast" );

        // no apoptosis until activation for homeostasis
        if( pCell.customData.get( "activated_immune_cell" ) < 0.5 )
        {
            phenotype.death.rates.set( apoptosis_index, 0d );//  [apoptosis_index] = 0.0;
        }
        else
        {
            phenotype.death.rates.set( apoptosis_index, pCD.phenotype.death.rates.get( apoptosis_index ) );//TODO: check
        }

        if( phenotype.death.dead )
        {
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;
            return;
        }
        pCell.phenotype.secretion.netExportRates[collagen_index] = ( ( ( 0.942 * TGF_beta ) / ( 0.174 + TGF_beta ) )
                * ( pCell.customData.get( "collagen_secretion_rate" ) ) * 2.52e-7 );


        for( int n = 0; n < pCell.getMicroenvironment().mesh.voxels.length; n++ )
        {
            double TGF_beta2 = pCell.getMicroenvironment().get( n )[antiinflammatory_cytokine_index]; //TODO: check

            if( TGF_beta2 > 0 )
            {
                pCell.customData.set( "activated_immune_cell", 1.0 );
                System.out.println( "Activated TGF" );
                return;
            }
        }
    }
}