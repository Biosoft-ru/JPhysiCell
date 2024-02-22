package ru.biosoft.physicell.sample_projects.interactions;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.BasicSignaling;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;

public class DifferentiatedPhenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // find my cell definition 
        CellDefinition pCD = CellDefinition.getCellDefinition( pCell.typeName );
        Microenvironment microenvironment = pCell.getMicroenvironment();
        // sample environment 
        int nR = microenvironment.findDensityIndex( "resource" );
        int nTox = microenvironment.findDensityIndex( "toxin" );
        int nDebris = microenvironment.findDensityIndex( "debris" );

        // if dead, release debris
        if( phenotype.death.dead == true )
        {
            phenotype.secretion.netExportRates[nDebris] = phenotype.volume.total;
            pCell.functions.updatePhenotype = null;
            return;
        }

        double[] samples = pCell.nearest_density_vector();
        double R = samples[nR];
        double toxin = samples[nTox];
        double signal = 0.0;
        double hill = 0.0;

        // pressure reduces proliferation 
        signal = pCell.state.simplePressure;
        double pressure_halfmax = pCD.custom_data.get( "cycling_pressure_halfmax" ); // 0.5 
        hill = BasicSignaling.Hill_response_function( signal, pressure_halfmax, 1.5 );
        double base_val = pCD.phenotype.cycle.data.getExitRate( 0 );

        phenotype.cycle.data.setExitRate( 0, ( 1 - hill ) * base_val );

        // resource reduces necrotic death 
        double max_val = 0.0028;
        int nNecrosis = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        // get same from bacteria
        double necrosis_saturation = pCD.custom_data.get( "necrosis_saturation_resource" ); // 0.075 
        double necrosis_threshold = pCD.custom_data.get( "necrosis_threshold_resource" ); // 0.15 

        phenotype.death.rates.set( nNecrosis,
                max_val * BasicSignaling.decreasing_linear_response_function( R, necrosis_saturation, necrosis_threshold ) );

        // toxin increases apoptotic death 
        int nApoptosis = phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );

        double toxicity_halfmax = pCD.custom_data.get( "toxicity_halfmax" ); // 0.2 
        double relative_max_tox_death = pCD.custom_data.get( "relative_max_toxicity" ); // 100 

        signal = toxin;
        base_val = pCD.phenotype.death.rates.get( nApoptosis );
        double max_response = base_val * relative_max_tox_death;
        hill = BasicSignaling.Hill_response_function( signal, toxicity_halfmax, 1.5 );
        // System.out.println( "tox: " + signal + " " + hill); 
        phenotype.death.rates.set( nApoptosis, base_val + ( max_response - base_val ) * hill );
    }

    @Override
    public String display()
    {
        return "Resource reduces necrosis. " + "Toxin increases apoptosis.";
    }
}