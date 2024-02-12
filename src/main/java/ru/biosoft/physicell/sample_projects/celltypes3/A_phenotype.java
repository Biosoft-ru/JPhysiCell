package ru.biosoft.physicell.sample_projects.celltypes3;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;

public class A_phenotype extends update_phenotype
{
    private Model model;

    public A_phenotype(Model model)
    {
        this.model = model;
    }

    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        CellDefinition pCD = CellDefinition.getCellDefinition( "A" );
        int nApoptosis = pCD.phenotype.death.find_death_model_index( "Apoptosis" );
        int nNecrosis = pCD.phenotype.death.find_death_model_index( "Necrosis" );

        if( phenotype.death.dead == true )
        {
            phenotype.secretion.setSecretionToZero();
            phenotype.secretion.setUptakeToZero();
            phenotype.motility.is_motile = false;
            pCell.functions.updatePhenotype = null;
            return;
        }

        // sample A, B, C, resource, and pressure 
        int nA = microenvironment.findDensityIndex( "signal A" );
        int nB = microenvironment.findDensityIndex( "signal B" );
        int nC = microenvironment.findDensityIndex( "signal C" );
        int nR = microenvironment.findDensityIndex( "resource" );

        double A = pCell.nearest_density_vector()[nA];
        double B = pCell.nearest_density_vector()[nB];
        double C = pCell.nearest_density_vector()[nC];
        double R = pCell.nearest_density_vector()[nR];
        double p = pCell.state.simplePressure;

        // necrotic death rate 
        double base_necrosis_rate = pCD.phenotype.death.rates.get( nNecrosis );
        double necrosis_threshold = model.getParameterDouble( "A_necrosis_threshold" );
        phenotype.death.rates.set( nNecrosis, 0.0 );

        if( R < necrosis_threshold )
        {
            phenotype.death.rates.set( nNecrosis, base_necrosis_rate * ( 1.0 - R / necrosis_threshold ) );
            //                phenotype.death.rates[nNecrosis] = base_necrosis_rate;
            //                phenotype.death.rates[nNecrosis] *= ( 1.0 - R / necrosis_threshold );
        }

        // cycle rate 
        double param0 = model.getParameterDouble( "A_base_cycle" ) * R;
        UpDownSignal sig = new UpDownSignal( model );
        sig.baseParameter = param0;
        sig.maxParameter = model.getParameterDouble( "A_max_cycle" );
        sig.addEffect( A, model.getParameter( "A_cycle_A" ) );// A 
        sig.addEffect( B, model.getParameter( "A_cycle_B" ) );// B
        sig.addEffect( C, model.getParameter( "A_cycle_C" ) );// C 
        phenotype.cycle.data.setTransitionRate( 0, 0, sig.computeEffect() );
        if( p > model.getParameterDouble( "A_cycle_pressure_threshold" ) )
        {
            phenotype.cycle.data.setTransitionRate( 0, 0, 0 );
        }

        // apoptotic rate 
        double base_death_rate = model.getParameterDouble( "A_base_death" );
        double max_death_rate = model.getParameterDouble( "A_max_death" );
        sig.reset();
        sig.baseParameter = base_death_rate;
        sig.maxParameter = max_death_rate;
        sig.addEffect( A, model.getParameter( "A_death_A" ) ); // A       
        sig.addEffect( B, model.getParameter( "A_death_B" ) ); // B         
        sig.addEffect( C, model.getParameter( "A_death_C" ) ); // C 
        sig.addEffect( C, model.getParameter( "A_death_R" ) ); // R 
        phenotype.death.rates.set( nApoptosis, sig.computeEffect() );
        if( p > model.getParameterDouble( "A_apoptosis_pressure_threshold" ) )
        {
            phenotype.death.rates.set( nApoptosis, 10.0 );
        }

        // speed 
        double base_speed = model.getParameterDouble( "A_base_speed" );
        double max_speed = model.getParameterDouble( "A_max_speed" );
        sig.reset();
        sig.baseParameter = base_speed;
        sig.maxParameter = max_speed;
        sig.addEffect( A, model.getParameter( "A_speed_A" ) );// A
        sig.addEffect( B, model.getParameter( "A_speed_B" ) );// B
        sig.addEffect( C, model.getParameter( "A_speed_C" ) );// C 
        sig.addEffect( C, model.getParameter( "A_speed_R" ) ); // R 
        phenotype.motility.migration_speed = sig.computeEffect();

        // secretion 
        double base_secretion = model.getParameterDouble( "A_base_secretion" );
        double max_secretion = model.getParameterDouble( "A_max_secretion" );
        sig.reset();
        sig.baseParameter = base_secretion;
        sig.maxParameter = max_secretion;
        sig.addEffect( A, model.getParameter( "A_signal_A" ) );// A            
        sig.addEffect( B, model.getParameter( "A_signal_B" ) );// B            
        sig.addEffect( C, model.getParameter( "A_signal_C" ) );// C            
        sig.addEffect( R, model.getParameter( "A_signal_R" ) ); // R
        phenotype.secretion.secretionRates[nA] = sig.computeEffect();
    }

    @Override
    public A_phenotype clone()
    {
        return new A_phenotype( model );
    }
}