package ru.biosoft.physicell.sample_projects.interactions;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.BasicSignaling;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;

public class StemPhenotype extends UpdatePhenotype
{
    double max_stem_diff;
    public StemPhenotype(Model model)
    {
        max_stem_diff = model.getParameterDouble( "max_stem_differentiation" );
    }
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

        // sample contacts 
        int stem_type = CellDefinition.getCellDefinition( "stem" ).type;
        int diff_type = CellDefinition.getCellDefinition( "differentiated" ).type;
        int bacteria_type = CellDefinition.getCellDefinition( "bacteria" ).type;

        int num_stem = 0;
        int num_differentiated = 0;
        int num_bacteria = 0;
        int num_dead = 0;
        //    for(
        //    int n = 0;n<pCell.state.neighbors.size();n++)
        //    {
        //        Cell pC = pCell.state.neighbors[n];
        for( Cell pC : pCell.state.neighbors )
        {
            if( pC.phenotype.death.dead == true )
            {
                num_dead++;
            }
            else
            {
                if( pC.type == stem_type )
                {
                    num_stem++;
                }
                if( pC.type == num_differentiated )
                {
                    num_differentiated++;
                }
                if( pC.type == bacteria_type )
                {
                    num_bacteria++;
                }
            }
        }

        // contact with a stem cell increases differentiation 
        //        double max_stem_diff = parameters.doubles( "max_stem_differentiation" ); // 0.0075 
        double stem_diff_halfmax = pCD.custom_data.get( "differentiation_contact_halfmax" ); // 0.1 

        double base_val = 0; // phenotype.cell_transformations.transformation_rates[diff_type]; 
        double max_val = max_stem_diff; // 0.0075;
        double signal = num_stem;
        double half_max = stem_diff_halfmax; // 0.1; 
        double hill = BasicSignaling.Hill_response_function( signal, half_max, 1.5 );
        phenotype.cellTransformations.transformationRates[diff_type] = base_val + ( max_val - base_val ) * hill;

        // contact with a differentiated cell reduces proliferation 
        // high rate of proliferation unless in contact with a differentiated cell 

        double stem_cycling_halfmax = pCD.custom_data.get( "cycling_contact_halfmax" ); // 0.1; 

        base_val = pCD.phenotype.cycle.data.getExitRate( 0 ); // 0.002; 
        max_val = 0.0;
        signal = num_differentiated;
        half_max = stem_cycling_halfmax; //  0.1; 
        hill = BasicSignaling.Hill_response_function( signal, half_max, 1.5 );
        phenotype.cycle.data.setExitRate( 0, base_val + ( max_val - base_val ) * hill );

        // resource reduces necrotic death 

        max_val = 0.0028;
        int nNecrosis = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
        double stem_saturation_necrosis = pCD.custom_data.get( "necrosis_saturation_resource" );
        double stem_threshold_necrosis = pCD.custom_data.get( "necrosis_threshold_resource" );

        phenotype.death.rates.set( nNecrosis,
                max_val * BasicSignaling.decreasing_linear_response_function( R, stem_saturation_necrosis, stem_threshold_necrosis ) );

        // toxin increases apoptotic death 
        int nApoptosis = phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );

        double toxicity_halfmax = pCD.custom_data.get( "toxicity_halfmax" ); // 0.4 
        double relative_max_toxicity = pCD.custom_data.get( "relative_max_toxicity" );

        signal = toxin;
        base_val = pCD.phenotype.death.rates.get( nApoptosis );
        max_val = base_val * relative_max_toxicity; // 100*base_val;

        hill = BasicSignaling.Hill_response_function( signal, toxicity_halfmax, 1.5 );
        phenotype.death.rates.set( nApoptosis, base_val + ( max_val - base_val ) * hill );
    }

    @Override
    public String display()
    {
        return "Contact with a stem cell increases differentiation." + " Contact with a differentiated cell reduces proliferation."
                + " Toxin increases apoptosis.";
    }
}