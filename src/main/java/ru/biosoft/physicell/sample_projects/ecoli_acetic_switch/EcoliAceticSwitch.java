package ru.biosoft.physicell.sample_projects.ecoli_acetic_switch;

import java.io.File;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;

import biouml.model.Diagram;
import biouml.plugins.fbc.ApacheModel;
import biouml.plugins.fbc.ApacheModelCreator;
import biouml.plugins.fbc.SbmlModelFBCReader2;
import biouml.plugins.sbml.SbmlModelFactory;
import biouml.plugins.sbml.SbmlModelReader_31;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.Secretion;
import ru.biosoft.physicell.core.Volume;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.core.standard.UpOrientation;
import ru.biosoft.physicell.fba.IntracellularFBA;

public class EcoliAceticSwitch extends Model
{
    CellDefinition bacteria_cell;
    IntracellularFBA fba;

    @Override
    public void init() throws Exception
    {
        super.init();
        setup_default_metabolic_model();
        createCellTypes();
        this.setupTissue();
        updater = new update_cell( this );
    }

    protected void doStep() throws Exception
    {
        //        m.simulateDiffusionDecay( diffusion_dt );
        ( (CellContainer)m.agentContainer ).updateAllCells( m, curTime, phenotype_dt, mechanics_dt, diffusion_dt );
        updateIntracellular();
        curTime += diffusion_dt;
        m.time = curTime;
    }

    void createBacteriaDefinition() throws Exception
    {
        bacteria_cell = StandardModels.getDefaultCellDefinition().clone( "bacteria cell", 1, m );
        bacteria_cell.phenotype.sync( m );
        bacteria_cell.phenotype.sync();
        CellDefinition.registerCellDefinition( bacteria_cell );
        double four_thirds_pi = 4.188790204786391;
        double cell_radius = getParameterDouble( "cell_radius" );
        double total = Math.pow( four_thirds_pi * cell_radius, 3 );
        double fluid_fraction = 0.75;
        double fluid = fluid_fraction * total;
        Volume volume = bacteria_cell.phenotype.volume;
        volume.total = total;
        volume.fluid = fluid;
        volume.solid = total - fluid;
        volume.nuclear = 0;
        volume.nuclear_fluid = 0;
        volume.nuclear_solid = 0;
        volume.cytoplasmic = total;
        volume.cytoplasmic_fluid = fluid_fraction * total;
        volume.cytoplasmic_solid = total - fluid;

        // Make sure we're ready for 2D
        bacteria_cell.functions.set_orientation = new UpOrientation();
        bacteria_cell.phenotype.geometry.polarity = 1.0;
        bacteria_cell.phenotype.motility.restrictTo2D = true;

        // use default proliferation and death
        int cycle_start_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
        int cycle_end_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
        int apoptosis_index = bacteria_cell.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );

        bacteria_cell.phenotype.death.rates.set( apoptosis_index, 0.0 );
        bacteria_cell.parameters.o2_proliferation_saturation = 38.0;
        bacteria_cell.parameters.o2_reference = 38.0;

        // set oxygen uptake and secretion to zero
        int oxygen_idx = m.findDensityIndex( "oxygen" ); // 0
        bacteria_cell.phenotype.secretion.secretionRates[oxygen_idx] = 0;
        bacteria_cell.phenotype.secretion.uptakeRates[oxygen_idx] = 0;
        bacteria_cell.phenotype.secretion.saturationDensities[oxygen_idx] = 0;

        int glucose_idx = m.findDensityIndex( "glucose" );
        bacteria_cell.phenotype.secretion.secretionRates[glucose_idx] = 0;
        bacteria_cell.phenotype.secretion.uptakeRates[glucose_idx] = 0;
        bacteria_cell.phenotype.secretion.saturationDensities[glucose_idx] = 0;

        int acetate_idx = m.findDensityIndex( "acetate" );
        bacteria_cell.phenotype.secretion.secretionRates[acetate_idx] = 0;
        bacteria_cell.phenotype.secretion.uptakeRates[acetate_idx] = 0;
        bacteria_cell.phenotype.secretion.saturationDensities[acetate_idx] = 0;

        // set the default cell type to no phenotype updates
        bacteria_cell.functions.updatePhenotype = null;
        bacteria_cell.functions.updateVolume = null;// = new update_cell( this );
    }

    void createCellTypes() throws Exception
    {
        CellDefinition cell_defaults = StandardModels.createDefaultCellDefinition();
        CellDefinition.registerCellDefinition( cell_defaults );
        cell_defaults.phenotype.secretion.sync( m );

        // turn the default cycle model to live, so it's easier to turn off proliferation
        cell_defaults.phenotype.cycle = StandardModels.live.clone();

        // Make sure we're ready for 2D
        cell_defaults.functions.set_orientation = new UpOrientation();
        cell_defaults.phenotype.geometry.polarity = 1.0;
        cell_defaults.phenotype.motility.restrictTo2D = true; // true; 

        // set to no motility for cancer cells 
        cell_defaults.phenotype.motility.isMotile = false;

        // turn the default cycle model to live,
        // so it's easier to turn off proliferation
        cell_defaults.name = "metabolic cell";
        cell_defaults.type = 0;

        createBacteriaDefinition();
    }

    //void setup_microenvironment( )
    //{
    //
    //	if( default_microenvironment_options.simulate_2D == false )
    //	{
    ////		std::cout << "WARNING: overriding from 3-D to 2-D" << std::endl;
    //		default_microenvironment_options.simulate_2D = true;
    //	}

    //	default_microenvironment_options.calculate_gradients = true;
    //	default_microenvironment_options.track_internalized_substrates_in_each_agent = false;

    //	initialize_microenvironment();
    //
    //	return;
    //}


    void setupTissue()
    {
        CellDefinition cd = CellDefinition.getCellDefinition( "bacteria cell" );
        // place a bacterial colony at the center 
        double cell_radius = cd.phenotype.geometry.radius;
        double cell_spacing = 0.95 * 2.0 * cell_radius;
        double colony_radius = 4;//getParameterDouble( "colony_radius" );

        Cell pCell = null;

        double x = 0.0;
        double x_outer = colony_radius;
        double y = 0.0;

        //        pCell = Cell.createCell( bacteria_cell, m, new double[] {0, 0, 0.0} );
        //        pCell.phenotype.intracellular = fba.clone();

        int n = 0;
        while( y < colony_radius )
        {
            x = 0.0;
            if( n % 2 == 1 )
            {
                x = 0.5 * cell_spacing;
            }
            x_outer = Math.sqrt( colony_radius * colony_radius - y * y );

            while( x < x_outer )
            {
                pCell = Cell.createCell( bacteria_cell, m, new double[] {x, y, 0.0} );
                pCell.phenotype.intracellular = fba.clone();

                if( Math.abs( y ) > 0.01 )
                {
                    pCell = Cell.createCell( bacteria_cell, m, new double[] {x, -y, 0.0} );
                    pCell.phenotype.intracellular = fba.clone();
                }

                if( Math.abs( x ) > 0.01 )
                {
                    pCell = Cell.createCell( bacteria_cell, m, new double[] { -x, y, 0.0} );
                    pCell.phenotype.intracellular = fba.clone();
                    if( Math.abs( y ) > 0.01 )
                    {
                        pCell = Cell.createCell( bacteria_cell, m, new double[] { -x, -y, 0.0} );
                        pCell.phenotype.intracellular = fba.clone();
                    }
                }
                x += cell_spacing;

            }

            y += cell_spacing * Math.sqrt( 3.0 ) / 2.0;
            n++;
        }
    }

    void setup_default_metabolic_model() throws Exception
    {
        fba = new IntracellularFBA();
        String path = "C:/Users/Damag/git/JPhysiCell/src/main/resources/ru/biosoft/physicell/sample_projects/ecoli_acetic_switch/config/Ecoli_core.xml";
        //        SbmlModelReader reader = SbmlModelFactory.readDiagram( path, false );
        Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse( new File( path ) );
        SbmlModelReader_31 reader = (SbmlModelReader_31)SbmlModelFactory.getReader( document );
        SbmlModelFBCReader2 packageReader = new SbmlModelFBCReader2();
        Diagram diagram = reader.read( document, "EcoliAceticSwitch", null, packageReader );
        ApacheModelCreator creator = new ApacheModelCreator();
        fba.model.fbcModel = (ApacheModel)creator.createModel( diagram );
        fba.parameterMapping = new HashMap<>();
        fba.parameterMapping.put( "glucose", "R_EX_glc__D_e" );
        fba.parameterMapping.put( "acetate", "R_EX_ac_e" );
        fba.parameterMapping.put( "oxygen", "R_EX_o2_e" );
        fba.parameterMapping.put( "growth_rate", "R_BIOMASS_Ecoli_core_w_GAM" );
    }

    void anuclear_volume_model(Cell pCell, Phenotype phenotype, double dt)
    {
        phenotype.volume.fluid += dt * phenotype.volume.fluid_change_rate
                * ( phenotype.volume.target_fluid_fraction * phenotype.volume.total - phenotype.volume.fluid );

        // if the fluid volume is negative, set to zero
        if( phenotype.volume.fluid < 0.0 )
        {
            phenotype.volume.fluid = 0.0;
        }


        phenotype.volume.cytoplasmic_solid += dt * phenotype.volume.cytoplasmic_biomass_change_rate
                * ( phenotype.volume.target_solid_cytoplasmic - phenotype.volume.cytoplasmic_solid );

        if( phenotype.volume.cytoplasmic_solid < 0.0 )
        {
            phenotype.volume.cytoplasmic_solid = 0.0;
        }

        phenotype.volume.solid = phenotype.volume.cytoplasmic_solid;
        phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + phenotype.volume.cytoplasmic_fluid;

        phenotype.volume.total = phenotype.volume.cytoplasmic;


        phenotype.volume.fluid_fraction = phenotype.volume.fluid / ( 1e-16 + phenotype.volume.total );

        phenotype.geometry.update( pCell, phenotype, dt );
    }

    //    void metabolic_cell_phenotype(Cell pCell, Phenotype phenotype, double dt)
    //    {
    //        // if cell is dead, don't bother with future phenotype changes.
    //        if( phenotype.death.dead == true )
    //        {
    //            pCell.functions.updatePhenotype = null;
    //            return;
    //        }
    //
    //        // update the transition rate according to growth rate?
    //        int cycle_start_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
    //        int cycle_end_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
    //
    //        //static int oncoprotein_i = pCell.custom_data.find_variable_index( "oncoprotein" );
    //        //phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell.custom_data[oncoprotein_i] ;
    //    }



    //    Color[] my_coloring_function(Cell pCell)
    //    {
    //        // start with flow cytometry coloring
    //
    //        List<String> output = false_cell_coloring_cytometry( pCell );
    //        output[0] = "red";
    //        output[1] = "red";
    //        output[2] = "red";
    //
    //        if( !pCell.phenotype.death.dead && pCell.type == 1 )
    //        {
    //            output[0] = "black";
    //            output[2] = "black";
    //        }
    //
    //        return output;
    //    }

    private update_cell updater;

    @Override
    public void updateIntracellular() throws Exception
    {
        for( Cell cell : m.getAgents( Cell.class ) )
            updater.execute( cell, cell.phenotype, curTime );
    }

    @Override
    public String getReportHeader()
    {
        return "ID\tX\tY\tZ\toxygen\tglucose\tacetate";
    }

    @Override
    public String getReport(Cell cell) throws Exception
    {
        Microenvironment m = cell.getMicroenvironment();
        int oxygen_idx = m.findDensityIndex( "oxygen" );
        int glucose_idx = m.findDensityIndex( "glucose" );
        int acetate_idx = m.findDensityIndex( "acetate" );

        Secretion secretion = cell.phenotype.secretion;
        double uptakeOxygen = secretion.uptakeRates[oxygen_idx];
        double uptakeGlucose = secretion.uptakeRates[glucose_idx];
        double uptakeAcetate = secretion.uptakeRates[acetate_idx];

        return "\n" + cell.ID + "\t" + cell.position[0] + "\t" + cell.position[1] + "\t" + cell.position[2] + "\t" + uptakeOxygen + "\t"
                + uptakeGlucose + "\t" + uptakeAcetate;
    }
}