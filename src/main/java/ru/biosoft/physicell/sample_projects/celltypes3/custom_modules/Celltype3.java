package ru.biosoft.physicell.sample_projects.celltypes3.custom_modules;

import java.awt.Color;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.ui.AgentVisualizer;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.xml.ModelReader;

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
public class Celltype3
{
    public static void init(Model model) throws Exception
    {
        SignalBehavior.setup_signal_behavior_dictionaries( model.getMicroenvironment() );
        create_cell_types( model );
        setup_tissue( model );
        for( Visualizer visualizer : model.getVisualizers() )
        {
            visualizer.setAgentVisualizer( new FluorescenceAgentVisualizer() );
        }
    }

    static void create_cell_types(Model model) throws Exception
    {
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );
        // set the random seed 
        //        SeedRandom( parameters.ints( "random_seed" ) );

        /* 
           Put any modifications to default cell definition here if you 
           want to have "inherited" by other cell types. 
           
           This is a good place to set default functions. 
        */

        //        initialize_default_cell_definition();
        CellDefinition cell_defaults = StandardModels.getDefaultCellDefinition();
        cell_defaults.functions.updateVolume = new StandardModels.standard_volume_update_function();
        cell_defaults.functions.updateVelocity = new StandardModels.standard_update_cell_velocity();

        cell_defaults.functions.update_migration_bias = null;
        cell_defaults.functions.updatePhenotype = null; // update_cell_and_death_parameters_O2_based; 
        cell_defaults.functions.custom_cell_rule = null;

        cell_defaults.functions.add_cell_basement_membrane_interactions = null;
        cell_defaults.functions.calculate_distance_to_membrane = null;



        CellDefinition.getCellDefinition( "A" ).functions.updatePhenotype = new A_phenotype( model );
        CellDefinition.getCellDefinition( "B" ).functions.updatePhenotype = new B_phenotype( model );
        CellDefinition.getCellDefinition( "C" ).functions.updatePhenotype = new C_phenotype( model );
    }

    static void setup_tissue(Model model)
    {
        Microenvironment microenvironment = model.getMicroenvironment();
        double Xmin = microenvironment.mesh.boundingBox[0];
        double Ymin = microenvironment.mesh.boundingBox[1];
        double Zmin = microenvironment.mesh.boundingBox[2];

        double Xmax = microenvironment.mesh.boundingBox[3];
        double Ymax = microenvironment.mesh.boundingBox[4];
        double Zmax = microenvironment.mesh.boundingBox[5];

        double max_radius = model.getParameterDouble( "max_distance_from_origin" );
        if( Xmax > max_radius )
        {
            Xmax = max_radius;
        }
        if( Xmin < -max_radius )
        {
            Xmin = -max_radius;
        }

        if( Ymax > max_radius )
        {
            Ymax = max_radius;
        }
        if( Ymin < -max_radius )
        {
            Ymin = -max_radius;
        }

        if( Zmax > max_radius )
        {
            Zmax = max_radius;
        }
        if( Zmin < -max_radius )
        {
            Zmin = -max_radius;
        }

        if( microenvironment.options.simulate_2D == true )
        {
            Zmin = 0.0;
            Zmax = 0.0;
        }

        double Xrange = Xmax - Xmin;
        double Yrange = Ymax - Ymin;
        double Zrange = Zmax - Zmin;

        //  Xmin += 0.25*Xrange; 
        //  Xrange *= 0.5;

        //  Ymin += 0.25*Yrange; 
        //  Yrange *= 0.5;

        //  Zmin += 0.25*Zrange; 
        //  Zrange *= 0.5; 
        // create some of each type of cell 

        // place A
        CellDefinition aCD = CellDefinition.getCellDefinition( "A" );
        CellDefinition bCD = CellDefinition.getCellDefinition( "B" );
        CellDefinition cCD = CellDefinition.getCellDefinition( "C" );

        for( int n = 0; n < model.getParameterInt( "number_of_A" ); n++ )
        {
            double[] position = {0, 0, 0};

            double r = max_radius + 1;
            while( r > max_radius )
            {
                position[0] = Xmin + PhysiCellUtilities.UniformRandom() * Xrange;
                position[1] = Ymin + PhysiCellUtilities.UniformRandom() * Yrange;
                position[2] = Zmin + PhysiCellUtilities.UniformRandom() * Zrange;

                r = VectorUtil.norm( position );
            }
            Cell pC = Cell.createCell( aCD, microenvironment, position );
            for( int k = 0; k < pC.phenotype.death.rates.size(); k++ )
            {
                pC.phenotype.death.rates.set( k, 0.0 );
            }
        }

        // place B
        for( int n = 0; n < model.getParameterInt( "number_of_B" ); n++ )
        {
            double[] position = {0, 0, 0};

            double r = max_radius + 1;
            while( r > max_radius )
            {
                position[0] = Xmin + PhysiCellUtilities.UniformRandom() * Xrange;
                position[1] = Ymin + PhysiCellUtilities.UniformRandom() * Yrange;
                position[2] = Zmin + PhysiCellUtilities.UniformRandom() * Zrange;

                r = VectorUtil.norm( position );
            }
            Cell pC = Cell.createCell( bCD, microenvironment, position );
            for( int k = 0; k < pC.phenotype.death.rates.size(); k++ )
            {
                pC.phenotype.death.rates.set( k, 0.0 );
            }
        }

        // place C
        for( int n = 0; n < model.getParameterInt( "number_of_C" ); n++ )
        {
            double[] position = {0, 0, 0};

            double r = max_radius + 1;
            while( r > max_radius )
            {
                position[0] = Xmin + PhysiCellUtilities.UniformRandom() * Xrange;
                position[1] = Ymin + PhysiCellUtilities.UniformRandom() * Yrange;
                position[2] = Zmin + PhysiCellUtilities.UniformRandom() * Zrange;

                r = VectorUtil.norm( position );
            }

            Cell pC = Cell.createCell( cCD, microenvironment, position );
            for( int k = 0; k < pC.phenotype.death.rates.size(); k++ )
            {
                pC.phenotype.death.rates.set( k, 0.0 );
            }
        }
    }

    public static class RegularAgentVisualizer extends AgentVisualizer
    {
        Model model;
        Color aColor;
        Color bColor;
        Color cColor;

        public RegularAgentVisualizer(Model model)
        {
            this.model = model;
            aColor = ModelReader.readColor( model.getParameter( "A_color" ) );
            bColor = ModelReader.readColor( model.getParameter( "B_color" ) );
            cColor = ModelReader.readColor( model.getParameter( "C_color" ) );
        }
        @Override
        public Color findColor(Cell pCell)
        {
            CellDefinition aCD = CellDefinition.getCellDefinition( "A" );
            CellDefinition bCD = CellDefinition.getCellDefinition( "B" );
            CellDefinition cCD = CellDefinition.getCellDefinition( "C" );

            int A_type = aCD.type;
            int B_type = bCD.type;
            int C_type = cCD.type;

            // start with flow cytometry coloring 
            Color output = Color.black;
            //            String[] output = {"black", "black", "black", "black"};

            // color live C 
            if( pCell.type == A_type )
            {
                output = aColor;//ModelReader.readColor( model.getParameter( "A_color" ) );
                //                output[0] = parameters.strings( "A_color" );
                //                output[2] = parameters.strings( "A_color" );
            }

            // color live B
            if( pCell.type == B_type )
            {
                output = bColor;//ModelReader.readColor( model.getParameter( "B_color" ) );
                //                output[0] = parameters.strings( "B_color" );
                //                output[2] = parameters.strings( "B_color" );
            }

            // color live C
            if( pCell.type == C_type )
            {
                output = cColor;//ModelReader.readColor( model.getParameter( "C_color" ) );
                //                output[0] = parameters.strings( "C_color" );
                //                output[2] = parameters.strings( "C_color" );
            }

            if( pCell.phenotype.death.dead == true )
            {
                // Necrotic - Brown
                if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                        || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                        || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
                {
                    output = new Color( 123, 63, 0 );//  "chocolate";
                }
                else
                {
                    output = Color.black;//"black";
                }
            }
            return output;
        }
    }

    public static class FluorescenceAgentVisualizer extends AgentVisualizer
    {
        @Override
        public Color findColor(Cell pCell)
        {
            Microenvironment microenvironment = pCell.getMicroenvironment();
            CellDefinition aCD = CellDefinition.getCellDefinition( "A" );
            CellDefinition bCD = CellDefinition.getCellDefinition( "B" );
            CellDefinition cCD = CellDefinition.getCellDefinition( "C" );
            int A_type = aCD.type;
            int B_type = bCD.type;
            int C_type = cCD.type;

            int nA = microenvironment.findDensityIndex( "signal A" );
            int nB = microenvironment.findDensityIndex( "signal B" );
            int nC = microenvironment.findDensityIndex( "signal C" );

            // start with flow cytometry coloring 
            Color result = Color.black;
            //        String[] result = {"black", "black", "black", "black"};

            double max_fluorescence = 1; // 
            double value = 0.0;

            // color live A
            if( pCell.type == A_type )
            {
                value = pCell.phenotype.secretion.secretionRates[nA] / ( 0.001 + aCD.phenotype.secretion.secretionRates[nA] );

                value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
                if( pCell.phenotype.death.dead == true )
                {
                    value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
                }
                result = new Color( 1.0f, 0.0f, 1.0f, (float)value );
                //            sprintf( color, "rgba(255,0,255,%f)", value );
            }

            // color live B
            if( pCell.type == B_type )
            {
                value = pCell.phenotype.secretion.secretionRates[nB] / ( 0.001 + bCD.phenotype.secretion.secretionRates[nB] );
                value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
                if( pCell.phenotype.death.dead == true )
                {
                    value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
                }
                result = new Color( 0.0f, 1.0f, 0.0f, (float)value );//new Color(0,255,0);
                //            result.
                //            sprintf( color, "rgba(0,255,0,%f)", value );
            }

            // color live C
            if( pCell.type == C_type )
            {
                value = pCell.phenotype.secretion.secretionRates[nC] / ( 0.001 + cCD.phenotype.secretion.secretionRates[nC] );
                value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
                if( pCell.phenotype.death.dead == true )
                {
                    value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
                }
                result = new Color( 0.0f, 1.0f, 1.0f, (float)value );
                //            sprintf( color, "rgba(0,255,255,%f)", value );
            }

            // Necrotic - black
            if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                    || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                    || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
            {
                result = new Color( 0.0f, 0.0f, 0.0f, (float)value );
                //            sprintf( color, "rgba(0,0,0,%f)", value );
            }

            //        String[] output = {color, "none", color, "none"};
            return result;
        }
    }

    private static class up_down_signal
    {
        double up;
        double down;

        boolean no_promoters;
        boolean no_inhibitors;

        double base_parameter;
        double max_parameter;
        private Model model;

        public up_down_signal(Model model)
        {
            this.model = model;
            up = 0.0;
            down = 0.0;

            base_parameter = 0.0;
            max_parameter = 1.0;

            no_promoters = true;
            no_inhibitors = true;
        }

        void add_effect(double factor, char factor_type)
        {
            // neutral signal 
            if( factor_type == 'N' || factor_type == 'n' )
            {
                return;
            }

            // promoter signal 
            if( factor_type == 'P' || factor_type == 'p' )
            {
                // up = sum of all (scaled) promoter signals 
                up += factor;
                no_promoters = false;
                return;
            }

            // inhibitor signal 
            if( factor_type == 'I' || factor_type == 'i' )
            {
                down += factor;
                no_inhibitors = false;
                return;
            }
        }

        void add_effect(double factor, String factor_type)
        {
            this.add_effect( factor, factor_type.charAt( 0 ) );
        }

        double compute_effect_hill()
        {
            double hill_power = model.getParameterDouble( "hill_power" );
            double half_max = model.getParameterDouble( "half_max" );
            double denom_constant = Math.pow( half_max, hill_power );

            double temp = Math.pow( up, hill_power );
            double UP = temp / ( denom_constant + temp );
            if( no_promoters )
            {
                UP = 0.0;
            }

            temp = Math.pow( down, hill_power );
            double DOWN = denom_constant / ( denom_constant + temp );
            if( no_inhibitors )
            {
                DOWN = 1.0;
            }

            // return UP * DOWN; 
            return ( base_parameter + ( max_parameter - base_parameter ) * UP ) * DOWN;
        }

        double compute_effect()
        {
            return this.compute_effect_hill();
        }

        void reset()
        {
            up = 0.0;
            down = 0.0;
            no_promoters = true;
            no_inhibitors = true;

            base_parameter = 0.0;
            max_parameter = 1.0;

            return;
        }


        public String toString()
        {
            StringBuilder sb = new StringBuilder();
            sb.append( "up    : " + up + " (no promoters : " + no_promoters + ")\n" );
            sb.append( "down  : " + down + " (no inhibiters: " + no_inhibitors + ")\n" );
            sb.append( "effect: " + compute_effect() );// << std::endl; 
            return sb.toString();
        }
    }

    public static class A_phenotype implements update_phenotype
    {
        private Model model;
        public A_phenotype(Model model)
        {
            this.model = model;
        }

        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            Microenvironment microenvironment = pCell.getMicroenvironment();
            // housekeeping 
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

            up_down_signal sig = new up_down_signal( model );
            sig.base_parameter = param0;
            sig.max_parameter = model.getParameterDouble( "A_max_cycle" );

            // A 
            sig.add_effect( A, model.getParameter( "A_cycle_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "A_cycle_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "A_cycle_C" ) );

            phenotype.cycle.data.setTransitionRate( 0, 0, sig.compute_effect() );
            //            .transition_rate( 0, 0 ) = sig.compute_effect();
            if( p > model.getParameterDouble( "A_cycle_pressure_threshold" ) )
            {
                phenotype.cycle.data.setTransitionRate( 0, 0, 0 );
                //                .transition_rate( 0, 0 ) = 0.0;
            }

            // apoptotic rate 

            double base_death_rate = model.getParameterDouble( "A_base_death" );
            double max_death_rate = model.getParameterDouble( "A_max_death" );
            sig.reset();
            sig.base_parameter = base_death_rate;
            sig.max_parameter = max_death_rate;

            // A 
            sig.add_effect( A, model.getParameter( "A_death_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "A_death_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "A_death_C" ) );
            // R 
            sig.add_effect( C, model.getParameter( "A_death_R" ) );

            phenotype.death.rates.set( nApoptosis, sig.compute_effect() );
            if( p > model.getParameterDouble( "A_apoptosis_pressure_threshold" ) )
            {
                phenotype.death.rates.set( nApoptosis, 10.0 );
            }

            // speed 
            double base_speed = model.getParameterDouble( "A_base_speed" );
            double max_speed = model.getParameterDouble( "A_max_speed" );
            sig.reset();
            sig.base_parameter = base_speed;
            sig.max_parameter = max_speed;

            // A 
            sig.add_effect( A, model.getParameter( "A_speed_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "A_speed_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "A_speed_C" ) );
            // R 
            sig.add_effect( C, model.getParameter( "A_speed_R" ) );

            phenotype.motility.migration_speed = sig.compute_effect();

            // secretion 
            double base_secretion = model.getParameterDouble( "A_base_secretion" );
            double max_secretion = model.getParameterDouble( "A_max_secretion" );
            sig.reset();
            sig.base_parameter = base_secretion;
            sig.max_parameter = max_secretion;
            // A 
            sig.add_effect( A, model.getParameter( "A_signal_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "A_signal_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "A_signal_C" ) );
            // R 
            sig.add_effect( R, model.getParameter( "A_signal_R" ) );

            phenotype.secretion.secretionRates[nA] = sig.compute_effect();
        }
    }

    public static class B_phenotype implements update_phenotype
    {
        private Model model;
        public B_phenotype(Model model)
        {
            this.model = model;
        }

        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // housekeeping 
            Microenvironment microenvironment = pCell.getMicroenvironment();
            CellDefinition pCD = CellDefinition.getCellDefinition( "B" );
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
            double necrosis_threshold = model.getParameterDouble( "B_necrosis_threshold" );
            phenotype.death.rates.set( nNecrosis, 0.0 );

            if( R < necrosis_threshold )
            {
                phenotype.death.rates.set( nNecrosis, base_necrosis_rate * ( 1.0 - R / necrosis_threshold ) );
                //                phenotype.death.rates[nNecrosis] = base_necrosis_rate;
                //                phenotype.death.rates[nNecrosis] *= ( 1.0 - R / necrosis_threshold );
            }

            // cycle rate 
            double param0 = model.getParameterDouble( "B_base_cycle" ) * R;

            up_down_signal sig = new up_down_signal( model );
            sig.base_parameter = param0;
            sig.max_parameter = model.getParameterDouble( "B_max_cycle" );

            // A 
            sig.add_effect( A, model.getParameter( "B_cycle_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "B_cycle_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "B_cycle_C" ) );

            //            phenotype.cycle.data.transition_rate( 0, 0 ) = sig.compute_effect();
            phenotype.cycle.data.setTransitionRate( 0, 0, sig.compute_effect() );
            if( p > model.getParameterDouble( "B_cycle_pressure_threshold" ) )
            {
                phenotype.cycle.data.setTransitionRate( 0, 0, 0 );
                //                phenotype.cycle.data.transition_rate( 0, 0 ) = 0.0;
            }

            // apoptotic rate 
            double base_death_rate = model.getParameterDouble( "B_base_death" );
            double max_death_rate = model.getParameterDouble( "B_max_death" );
            sig.reset();
            sig.base_parameter = base_death_rate;
            sig.max_parameter = max_death_rate;

            // A 
            sig.add_effect( A, model.getParameter( "B_death_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "B_death_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "B_death_C" ) );
            // R 
            sig.add_effect( C, model.getParameter( "B_death_R" ) );

            phenotype.death.rates.set( nApoptosis, sig.compute_effect() );
            if( p > model.getParameterDouble( "A_apoptosis_pressure_threshold" ) )
            {
                phenotype.death.rates.set( nApoptosis, 10.0 );
            }

            // speed 
            double base_speed = model.getParameterDouble( "B_base_speed" );
            double max_speed = model.getParameterDouble( "B_max_speed" );
            sig.reset();
            sig.base_parameter = base_speed;
            sig.max_parameter = max_speed;
            // A 
            sig.add_effect( A, model.getParameter( "B_speed_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "B_speed_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "B_speed_C" ) );
            // R 
            sig.add_effect( C, model.getParameter( "B_speed_R" ) );

            phenotype.motility.migration_speed = sig.compute_effect();

            // secretion 

            double base_secretion = model.getParameterDouble( "B_base_secretion" );
            double max_secretion = model.getParameterDouble( "B_max_secretion" );
            sig.reset();
            sig.base_parameter = base_secretion;
            sig.max_parameter = max_secretion;

            // A 
            sig.add_effect( A, model.getParameter( "B_signal_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "B_signal_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "B_signal_C" ) );
            // R 
            sig.add_effect( R, model.getParameter( "B_signal_R" ) );

            phenotype.secretion.secretionRates[nB] = sig.compute_effect();

        }
    }

    public static class C_phenotype implements update_phenotype
    {
        private Model model;
        public C_phenotype(Model model)
        {
            this.model = model;
        }

        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            Microenvironment microenvironment = pCell.getMicroenvironment();
            // housekeeping 
            CellDefinition pCD = CellDefinition.getCellDefinition( "C" );
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
            double param0 = model.getParameterDouble( "C_base_cycle" ) * R;

            up_down_signal sig = new up_down_signal( model );
            sig.base_parameter = param0;
            sig.max_parameter = model.getParameterDouble( "C_max_cycle" );

            // A 
            sig.add_effect( A, model.getParameter( "C_cycle_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "C_cycle_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "C_cycle_C" ) );

            //            phenotype.cycle.data.transition_rate( 0, 0 ) = sig.compute_effect();
            phenotype.cycle.data.setTransitionRate( 0, 0, sig.compute_effect() );
            if( p > model.getParameterDouble( "C_cycle_pressure_threshold" ) )
            {
                phenotype.cycle.data.setTransitionRate( 0, 0, 0 );
                //                phenotype.cycle.data.transition_rate( 0, 0 ) = 0.0;
            }

            // apoptotic rate 
            double base_death_rate = model.getParameterDouble( "C_base_death" );
            double max_death_rate = model.getParameterDouble( "C_max_death" );
            sig.reset();
            sig.base_parameter = base_death_rate;
            sig.max_parameter = max_death_rate;

            // A 
            sig.add_effect( A, model.getParameter( "C_death_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "C_death_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "C_death_C" ) );
            // R 
            sig.add_effect( C, model.getParameter( "C_death_R" ) );

            phenotype.death.rates.set( nApoptosis, sig.compute_effect() );
            if( p > model.getParameterDouble( "C_apoptosis_pressure_threshold" ) )
            {
                phenotype.death.rates.set( nApoptosis, 10.0 );
            }

            // speed 
            double base_speed = model.getParameterDouble( "C_base_speed" );
            double max_speed = model.getParameterDouble( "C_max_speed" );
            sig.reset();
            sig.base_parameter = base_speed;
            sig.max_parameter = max_speed;

            // A 
            sig.add_effect( A, model.getParameter( "C_speed_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "C_speed_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "C_speed_C" ) );
            // R 
            sig.add_effect( C, model.getParameter( "C_speed_R" ) );

            phenotype.motility.migration_speed = sig.compute_effect();

            // secretion 
            double base_secretion = model.getParameterDouble( "C_base_secretion" );
            double max_secretion = model.getParameterDouble( "C_max_secretion" );
            sig.reset();
            sig.base_parameter = base_secretion;
            sig.max_parameter = max_secretion;

            // A 
            sig.add_effect( A, model.getParameter( "C_signal_A" ) );
            // B
            sig.add_effect( B, model.getParameter( "C_signal_B" ) );
            // C 
            sig.add_effect( C, model.getParameter( "C_signal_C" ) );
            // R 
            sig.add_effect( R, model.getParameter( "C_signal_R" ) );

            phenotype.secretion.secretionRates[nC] = sig.compute_effect();
        }
    }
}