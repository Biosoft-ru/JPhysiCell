package ru.biosoft.physicell.sample_projects.interactions;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.BasicSignaling;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;

public class BacteriaPhenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( pCell.typeName );

        Microenvironment m = pCell.getMicroenvironment();
        // sample resource, quorum, and toxin 
        int nR = m.findDensityIndex( "resource" );
        int nDebris = m.findDensityIndex( "debris" );
        int nQuorum = m.findDensityIndex( "quorum" );
        int nToxin = m.findDensityIndex( "toxin" );

        // if dead: stop exporting quorum factor. 
        // also, replace phenotype function 
        if( phenotype.death.dead == true )
        {
            phenotype.secretion.netExportRates[nQuorum] = 0;
            phenotype.secretion.netExportRates[nToxin] = 0;

            phenotype.secretion.netExportRates[nDebris] = phenotype.volume.total;

            pCell.functions.updatePhenotype = null;
            return;
        }

        double[] samples = pCell.nearest_density_vector();
        double R = samples[nR];
        double Q = samples[nQuorum];
        double Tox = samples[nToxin];

        // resource increases cycle entry 
        double base_val = pCD.phenotype.cycle.data.getExitRate( 0 );
        double max_val = base_val * 10.0;
        double min_cycle_resource = pCD.custom_data.get( "cycling_entry_threshold_resource" ); // 0.15 
        phenotype.cycle.data.setExitRate( 0, max_val * BasicSignaling.linear_response_function( R, min_cycle_resource, 1 ) );

        // resource decreses necrosis
        max_val = 0.0028;
        int nNecrosis = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
        double saturation_necrosis_resource = pCD.custom_data.get( "necrosis_saturation_resource" ); //0.075
        double threshold_necrosis_resource = pCD.custom_data.get( "necrosis_threshold_resource" ); // 0.15
        phenotype.death.rates.set( nNecrosis, max_val
                * BasicSignaling.decreasing_linear_response_function( R, saturation_necrosis_resource, threshold_necrosis_resource ) );

        // resource decreases motile speed  
        double signal = R;
        base_val = pCD.phenotype.motility.migrationSpeed;
        double max_response = 0.0;
        double motility_resource_halfmax = pCD.custom_data.get( "migration_speed_halfmax" ); // 0.25 // parameters.doubles("bacteria_motility_resource_halfmax");
        double hill = BasicSignaling.Hill_response_function( signal, motility_resource_halfmax, 1.5 );
        phenotype.motility.migrationSpeed = base_val + ( max_response - base_val ) * hill;

        // quorum and resource increases motility bias 
        signal = Q + R;
        base_val = pCD.phenotype.motility.migrationSpeed;
        max_response = 1.0;
        double bias_halfmax = pCD.custom_data.get( "migration_bias_halfmax" );
        // 0.5 //  parameters.doubles("bacteria_migration_bias_halfmax");
        hill = BasicSignaling.Hill_response_function( signal, bias_halfmax, 1.5 );
        phenotype.motility.migrationBias = base_val + ( max_response - base_val ) * hill;

        // damage increases death 
        int nApoptosis = phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );

        signal = pCell.state.damage;
        base_val = pCD.phenotype.death.rates.get( nApoptosis );

        double damage_halfmax = pCD.custom_data.get( "damage_halfmax" );
        double relative_max_damage_death = pCD.custom_data.get( "relative_max_damage_death" );
        max_response = base_val * relative_max_damage_death;

        // 36 // parameters.doubles("bacteria_damage_halfmax");
        hill = BasicSignaling.Hill_response_function( signal, damage_halfmax, 1.5 );
        phenotype.death.rates.set( nApoptosis, base_val + ( max_response - base_val ) * hill );
    }

    @Override
    public String display()
    {
        return "Resource decrease necrosis and motility, also stimulate division." + " Resource AND quorum increase motility."
                + " Damage stimulates apoptosis";
    }
}