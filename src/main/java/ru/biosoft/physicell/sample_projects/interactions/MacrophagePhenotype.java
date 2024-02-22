package ru.biosoft.physicell.sample_projects.interactions;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.BasicSignaling;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;

/* https://www.karger.com/Article/Fulltext/494069 */
public class MacrophagePhenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {

        // find my cell definition 
        CellDefinition pCD = CellDefinition.getCellDefinition( pCell.typeName );
        Microenvironment microenvironment = pCell.getMicroenvironment();

        // sample environment 
        int nPIF = microenvironment.findDensityIndex( "pro-inflammatory" );
        int nDebris = microenvironment.findDensityIndex( "debris" );
        int nQ = microenvironment.findDensityIndex( "quorum" );

        // if dead, release debris
        if( phenotype.death.dead == true )
        {
            phenotype.secretion.netExportRates[nDebris] = phenotype.volume.total;
            pCell.functions.updatePhenotype = null;
            return;
        }

        double[] samples = pCell.nearest_density_vector();
        double PIF = samples[nPIF];
        double debris = samples[nDebris];
        double Q = samples[nQ];

        // sample contacts 
        int bacteria_type = CellDefinition.getCellDefinition( "bacteria" ).type;

        int num_bacteria = 0;
        int num_dead = 0;

        for( Cell pC : pCell.state.neighbors )
        {
            if( pC.phenotype.death.dead )
                num_dead++;
            else
            {
                if( pC.type == bacteria_type )
                    num_bacteria++;
            }
        }

        // contact with dead cells or bacteria, or debris increases secretion of pro-inflammatory 
        double secretion_dead_sensitivity = 1;
        double secretion_bacteria_sensitivity = 1;
        double secretion_debris_sensitivity = 2;
        double secretion_quorum_sensitivity = 5;

        double base_val = pCD.phenotype.secretion.secretionRates[nPIF];
        double max_response = 10; // phenotype.volume.total; 
        double signal = secretion_dead_sensitivity * num_dead + secretion_bacteria_sensitivity * num_bacteria
                + secretion_debris_sensitivity * debris + secretion_quorum_sensitivity * Q;
        double half_max = pCD.custom_data.get( "secretion_halfmax" ); // 0.5; // 0.5; 
        double hill = BasicSignaling.Hill_response_function( signal, half_max, 1.5 );
        phenotype.secretion.secretionRates[nPIF] = base_val + ( max_response - base_val ) * hill;

        // chemotaxis bias increases with debris or quorum factor 
        double bias_debris_sensitivity = 0.1;
        double bias_quorum_sensitivity = 1;

        base_val = pCD.phenotype.motility.migrationBias;
        max_response = 0.75;
        signal = bias_debris_sensitivity * debris + bias_quorum_sensitivity * Q; // + 10 * PIF; 
        half_max = pCD.custom_data.get( "migration_bias_halfmax" ); // 0.01 // 0.005 //0.1 // 0.05
        hill = BasicSignaling.Hill_response_function( signal, half_max, 1.5 );
        phenotype.motility.migrationBias = base_val + ( max_response - base_val ) * hill;

        // migration speed slows down in the presence of debris or quorum factor 
        base_val = pCD.phenotype.motility.migrationSpeed;
        max_response = 0.1 * base_val;
        signal = bias_debris_sensitivity * debris + bias_quorum_sensitivity * Q; // + 10 * PIF; 
        half_max = pCD.custom_data.get( "migration_speed_halfmax" ); // 0.1 // 0.05 
        hill = BasicSignaling.Hill_response_function( signal, half_max, 1.5 );
        phenotype.motility.migrationSpeed = base_val + ( max_response - base_val ) * hill;
    }

    @Override
    public String display()
    {
        return "Contact with dead cells, bacteria or debris increases secretion of pro-inflammatory. "
                + "Debris and quorum increase chemotaxis bias. " + "Migration speed slows down in the presence of debris or quorum factor.";
    }
}