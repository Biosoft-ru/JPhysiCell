package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import java.awt.Color;
import java.util.Set;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.contact_function;
import ru.biosoft.physicell.core.CellFunctions.custom_cell_rule;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Model.Event;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.ui.AgentVisualizer;
import ru.biosoft.physicell.ui.Visualizer;

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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

public class CancerBiorobots
{

    public static void init(Model model) throws Exception
    {
        SignalBehavior.setup_signal_behavior_dictionaries( model.getMicroenvironment() );
        create_cell_types( model );
        setup_tissue( model );
        model.addEvent( new TherapyEvent( model.getParameterDouble( "therapy_activation_time" ) ) );
        for( Visualizer visualizer : model.getVisualizers() )
        {
            visualizer.setAgentVisualizer( new CancerBiorobotsVisualizer() );
        }
    }

    private static void create_cell_types(Model m)
    {
        PhysiCellUtilities.setSeed( m.getParameterInt( "random_seed" ) );

        //cancer cell
        CellDefinition pCD = CellDefinition.getCellDefinition( "cancer cell" );
        pCD.functions.updatePhenotype = new tumor_cell_phenotype_with_therapy();
        pCD.parameters.o2_proliferation_saturation = 38.0;
        pCD.parameters.o2_reference = 38.0;

        // cargo cells
        pCD = CellDefinition.getCellDefinition( "cargo cell" );
        // figure out mechanics parameters
        pCD.phenotype.mechanics.relative_maximum_attachment_distance = pCD.custom_data.get( "max_attachment_distance" )
                / pCD.phenotype.geometry.radius;
        pCD.phenotype.mechanics.relative_detachment_distance = pCD.custom_data.get( "max_elastic_displacement" )
                / pCD.phenotype.geometry.radius;
        pCD.phenotype.mechanics.attachment_elastic_constant = pCD.custom_data.get( "elastic_coefficient" );

        // set functions
        pCD.functions.updatePhenotype = new cargo_cell_phenotype_rule();
        pCD.functions.custom_cell_rule = new cargo_cell_rule();
        pCD.functions.contact_function = new biorobots_contact_function();
        pCD.functions.update_migration_bias = null;

        // worker cells
        pCD = CellDefinition.getCellDefinition( "worker cell" );
        pCD.phenotype.mechanics.relative_maximum_attachment_distance = pCD.custom_data.get( "max_attachment_distance" )
                / pCD.phenotype.geometry.radius;

        pCD.phenotype.mechanics.relative_detachment_distance = pCD.custom_data.get( "max_elastic_displacement" )
                / pCD.phenotype.geometry.radius;
        pCD.phenotype.mechanics.attachment_elastic_constant = pCD.custom_data.get( "elastic_coefficient" );
        pCD.functions.updatePhenotype = null; // worker_cell_rule;
        pCD.functions.custom_cell_rule = new worker_cell_rule();
        pCD.functions.contact_function = new biorobots_contact_function();

        /*
         * This builds the map of cell definitions and summarizes the setup.
         */

        //	display_CellDefinitions(std::cout);
    }

    private static void setup_tissue(Model model) throws Exception
    {
        Microenvironment microenvironment = model.getMicroenvironment();
        double Xmin = microenvironment.mesh.boundingBox[0];
        double Ymin = microenvironment.mesh.boundingBox[1];
        double Zmin = microenvironment.mesh.boundingBox[2];
        double Xmax = microenvironment.mesh.boundingBox[3];
        double Ymax = microenvironment.mesh.boundingBox[4];
        double Zmax = microenvironment.mesh.boundingBox[5];

        if( microenvironment.options.simulate_2D == true )
        {
            Zmin = 0.0;
            Zmax = 0.0;
        }
        // create some of each type of cell

        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            System.out.println( "Placing cells of type " + cd.name + " ... " );
            for( int n = 0; n < model.getParameterInt( "number_of_cells" ); n++ )
            {
                double[] position = {0, 0, 0};
                position[0] = PhysiCellUtilities.UniformRandom( Xmin, Xmax );// * Xrange;
                position[1] = PhysiCellUtilities.UniformRandom( Ymin, Ymax );// + PhysiCellUtilities.UniformRandom() * Yrange;
                position[2] = PhysiCellUtilities.UniformRandom( Zmin, Zmax );// + PhysiCellUtilities.UniformRandom() * Zrange;
                Cell.createCell( cd, microenvironment, position );
            }
        }
        //	std::cout << std::endl; 

        // custom placement, place a cluster of tumor cells at the center
        CellDefinition defaults = StandardModels.getDefaultCellDefinition();
        double cell_radius = defaults.phenotype.geometry.radius;
        double cell_spacing = 0.95 * 2.0 * cell_radius;

        double tumor_radius = model.getParameterDouble( "tumor_radius" ); // 200.0;

        Cell pCell = null;
        CellDefinition pCD_cancer = CellDefinition.getCellDefinition( "cancer cell" );

        double x = 0.0;
        double x_outer = tumor_radius;
        double y = 0.0;

        int n = 0;
        while( y < tumor_radius )
        {
            x = 0.0;
            if( n % 2 == 1 )
            {
                x = 0.5 * cell_spacing;
            }
            x_outer = Math.sqrt( tumor_radius * tumor_radius - y * y );

            while( x < x_outer )
            {
                pCell = Cell.createCell( pCD_cancer, microenvironment, new double[] {x, y, 0.0} ); // tumor cell

                if( Math.abs( y ) > 0.01 )
                {
                    pCell = Cell.createCell( pCD_cancer, microenvironment, new double[] {x, -y, 0.0} ); // tumor cell			
                }

                if( Math.abs( x ) > 0.01 )
                {
                    Cell.createCell( pCD_cancer, microenvironment, new double[] { -x, y, 0.0} );

                    if( Math.abs( y ) > 0.01 )
                    {
                        Cell.createCell( pCD_cancer, microenvironment, new double[] { -x, -y, 0.0} );
                    }
                }
                x += cell_spacing;
            }

            y += cell_spacing * Math.sqrt( 3.0 ) / 2.0;
            n++;
        }
    }

    public static class CancerBiorobotsVisualizer extends AgentVisualizer
    {
        @Override
        public Color findBorderColor(Cell cell)
        {
            return Color.black;
        }

        @Override
        public Color findColor(Cell cell)
        {
            Color c = Color.white;
            double damage = SignalBehavior.get_single_signal( cell, "damage");
            double max_damage = 1.0 * SignalBehavior.get_single_signal(cell, "custom:damage_rate")
                    / ( 1e-16 + SignalBehavior.get_single_signal( cell, "custom:repair_rate" ) );

            CellDefinition pCD_cargo = CellDefinition.getCellDefinition( "cargo cell" );
            CellDefinition pCD_cancer = CellDefinition.getCellDefinition( "cancer cell" );
            CellDefinition pCD_worker = CellDefinition.getCellDefinition( "worker cell" );

            //   cargo cell 
            if( cell.type == pCD_cargo.type )
            {
                return Color.BLUE;
                //                output[0] = "blue";
                //                output[1] = "blue";
                //                output[2] = "blue";
                //                output[3] = "none"; // no nuclear outline color 
                //                return output;
            }
                    //  
                    // worker cell 
                    if( cell.type == pCD_worker.type )
                    {
                        return Color.red;
                    }
                    //      output[0] = "red";
                    //      output[1] = "red";
                    //      output[2] = "red"; 
                    //      output[3] = "none"; // no nuclear outline color 
                    //      return output;
                    //  }
                    //  
                    // apoptotic tumor - cyan 
                    if( SignalBehavior.get_single_signal( cell, "apoptotic" ) > 0.5 ) // Apoptotic - cyan
                    {
                        //                          output[0] = "cyan";
                        //                          output[2] = "darkcyan"; 
                        return Color.cyan;
                    }
                    //  
                    // Necrotic tumor - Brown
                    if( SignalBehavior.get_single_signal( cell, "necrotic" ) > 0.5 )
                    {
//                        output[0] = "rgb(250,138,38)";
//                        output[2] = "rgb(139,69,19)";
                        return new Color(250,138,18);
                    }
                    //  
                      // live tumor -- shade by level of damage 
                      // if live: color by damage 
                      if( SignalBehavior.get_single_signal( cell, "dead") < 0.5 )
                      {
                          int damage_int = (int) Math.round( damage * 255.0 / max_damage ); 
                          return new Color( damage_int, 255 - damage_int, damage_int );
                      }
            return c;
        }
    }

    public static void introduce_biorobots(Model model) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        // idea: we'll "inject" them in a little column
        double worker_fraction = model.getParameterDouble( "worker_fraction" ); // 0.10; /* param */
        int number_of_injected_cells = model.getParameterInt( "number_of_injected_cells" ); // 500; /* param */

        // make these vary with domain size
        double left_coordinate = m.options.X_range[1] - 150.0; // 600.0;
        double right_cooridnate = m.options.X_range[1] - 50.0; // 700.0;

        double bottom_coordinate = m.options.Y_range[0] + 50.0; // -700;
        double top_coordinate = m.options.Y_range[1] - 50.0; // 700;

        CellDefinition pCD_worker = CellDefinition.getCellDefinition( "worker cell" );
        CellDefinition pCD_cargo = CellDefinition.getCellDefinition( "cargo cell" );

        for( int i = 0; i < number_of_injected_cells; i++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( left_coordinate, right_cooridnate );//left_coordinate + ( right_cooridnate - left_coordinate ) * PhysiCellUtilities.UniformRandom();
            position[1] = PhysiCellUtilities.UniformRandom( bottom_coordinate, top_coordinate );//bottom_coordinate + ( top_coordinate - bottom_coordinate ) * PhysiCellUtilities.UniformRandom();

            if( PhysiCellUtilities.UniformRandom() <= worker_fraction )
            {
                Cell.createCell( pCD_worker, m, position );
            }
            else
            {
                Cell.createCell( pCD_cargo, m, position );
            }
        }
    }

    public static class cargo_cell_rule implements custom_cell_rule
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
        {
            if( SignalBehavior.get_single_signal( pCell, "dead" ) > 0.5 )
            {
                // the cell death functions don't automatically turn off custom functions, since those are part of mechanics.
                // Let's just fully disable now.
                pCell.functions.custom_cell_rule = null;
                return;
            }

            // if I'm docked
            if( pCell.state.numberAttachedCells() > 0 )
            {
                SignalBehavior.setSingleBehavior( pCell, "migration speed", 0.0 );
                return;
            }
        }
    }

    public static class cargo_cell_phenotype_rule extends update_phenotype
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
        {
            // if dettached and receptor on, secrete signal
            // if dettached and receptor off, secrete chemo
            double receptor = SignalBehavior.get_single_signal( pCell, "custom:receptor" );

            if( pCell.state.numberAttachedCells() == 0 )
            {
                if( receptor > 0.1 )
                {
                    SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 10 );
                    SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 0 );
                }
                else
                {
                    SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 0 );
                    SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 10 );
                }
                return;
            }

            // if you reach this point of the code, the cell is attached
            // if attached and oxygen high, secrete nothing, receptor off
            // if attached and oxygen low, dettach, start secreting chemo, receptor off
            double o2 = SignalBehavior.get_single_signal( pCell, "oxygen" );
            double o2_drop = SignalBehavior.get_single_signal( pCell, "custom:cargo_release_o2_threshold" );

            if( o2 > o2_drop )
            {
                SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 0 );
                SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 0 );
                SignalBehavior.setSingleBehavior( pCell, "custom:receptor", 0 );
            }
            else
            {
                SignalBehavior.setSingleBehavior( pCell, "chemoattractant secretion", 0 );
                SignalBehavior.setSingleBehavior( pCell, "therapeutic secretion", 10 );
                SignalBehavior.setSingleBehavior( pCell, "custom:receptor", 0 );
                pCell.remove_all_attached_cells();
            }
        }
    }

    public static class biorobots_contact_function implements contact_function
    {
        public void execute(Cell pActingOn, Phenotype pao, Cell pAttachedTo, Phenotype pat, double dt)
        {
            double[] displacement = VectorUtil.newDiff( pAttachedTo.position, pActingOn.position );
            double max_elastic_displacement = pao.geometry.radius * pao.mechanics.relative_detachment_distance;
            double max_displacement_squared = max_elastic_displacement * max_elastic_displacement;
            // detach cells if too far apart
            if( VectorUtil.norm_squared( displacement ) > max_displacement_squared )
            {
                Cell.detach_cells( pActingOn, pAttachedTo );
                return;
            }
            VectorUtil.axpy( ( pActingOn.velocity ), pao.mechanics.attachment_elastic_constant, displacement );
        }
    }

    public static class tumor_cell_phenotype_with_therapy extends update_phenotype
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
        {
            double damage = SignalBehavior.get_single_signal( pCell, "damage" );

            double damage_rate = SignalBehavior.get_single_signal( pCell, "custom:damage_rate" );
            double repair_rate = SignalBehavior.get_single_signal( pCell, "custom:repair_rate" );
            double drug_death_rate = SignalBehavior.get_single_signal( pCell, "custom:drug_death_rate" );

            double drug = SignalBehavior.get_single_signal( pCell, "therapeutic" );

            double max_damage = 1.0 * damage_rate / ( 1e-16 + repair_rate );

            // if I'm dead, don't bother. disable my phenotype rule
            if( SignalBehavior.get_single_signal( pCell, "dead" ) > 0.5 )
            {
                pCell.functions.updatePhenotype = null;
                return;
            }

            // first, vary the cell birth and death rates with oxygenation

            // std::cout << get_single_behavior( pCell , "cycle entry") << " vs ";
            new StandardModels.update_cell_and_death_parameters_O2_based().execute( pCell, phenotype, dt );
            // std::cout << get_single_behavior( pCell , "cycle entry") << std::endl;

            // the update the cell damage

            // dD/dt = alpha*c - beta-D by implicit scheme

            double temp = drug;

            // reuse temp as much as possible to reduce memory allocations etc.
            temp *= dt;
            temp *= damage_rate;

            damage += temp; // d_prev + dt*chemo*damage_rate

            temp = repair_rate;
            temp *= dt;
            temp += 1.0;
            damage /= temp; // (d_prev + dt*chemo*damage_rate)/(1 + dt*repair_rate)

            // then, see if the cell undergoes death from the therapy
            temp = dt;
            temp *= damage;
            temp *= drug_death_rate;
            temp /= max_damage; // dt*(damage/max_damage)*death_rate

            // make sure we write the damage (not current a behavior)
            pCell.state.damage = damage;
            if( damage > 0 )
            {
                System.out.println( damage );
            }
            if( PhysiCellUtilities.UniformRandom() <= temp )
            {
                // pCell.start_death( apoptosis_model_index );
                SignalBehavior.setSingleBehavior( pCell, "apoptosis", 9e99 );
                pCell.functions.updatePhenotype = null;
                pCell.functions.custom_cell_rule = null;
            }
        }
    }

    public static class worker_cell_rule implements custom_cell_rule
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
        {
            // if I am dead, don't bother

            if( SignalBehavior.get_single_signal( pCell, "dead" ) > 0.5 )
            {
                // the cell death functions don't automatically turn off custom functions,
                // since those are part of mechanics.

                // Let's just fully disable now.
                pCell.functions.custom_cell_rule = null;
                return;
            }

            // am I searching for cargo? if so, see if I've found it
            if( pCell.state.numberAttachedCells() == 0 )
            {
                Set<Cell> nearby = pCell.cells_in_my_container();
                boolean attached = false; // want to limit to one attachment
                for( Cell nearbyCell : nearby )
                {
                    if( nearbyCell == pCell )
                        continue;
                    // if it is expressing the receptor, dock with it
                    if( SignalBehavior.get_single_signal( nearbyCell, "custom:receptor" ) > 0.5 && attached == false )
                    {
                        Cell.attach_cells( pCell, nearbyCell );
                        // nearby[i].custom_data["receptor"] = 0.0; // put into cargo cell rule instead?
                        // nearby[i].phenotype.secretion.set_all_secretion_to_zero(); // put into cargo rule instead?
                        attached = true;
                        break;
                    }
                }
            }

            // from prior motility function
            double o2 = SignalBehavior.get_single_signal( pCell, "oxygen" );
            double chemoattractant = SignalBehavior.get_single_signal( pCell, "chemoattractant" );
            double detection_threshold = SignalBehavior.get_single_signal( pCell, "custom:motility_shutdown_detection_threshold" );

            // if attached, biased motility towards director chemoattractant
            // otherwise, biased motility towards cargo chemoattractant
            double attached_worker_migration_bias = SignalBehavior.get_single_signal( pCell, "custom:attached_worker_migration_bias" );
            double unattached_worker_migration_bias = SignalBehavior.get_single_signal( pCell, "custom:unattached_worker_migration_bias" );

            if( pCell.state.numberAttachedCells() > 0 )
            {
                SignalBehavior.setSingleBehavior( pCell, "migration bias", attached_worker_migration_bias );
                SignalBehavior.setSingleBehavior( pCell, "chemotactic response to oxygen", -1 );
                SignalBehavior.setSingleBehavior( pCell, "chemotactic response to chemoattractant", 0 );
            }
            else
            {
                // if there is no detectable signal, shut down motility (permanently)
                if( chemoattractant < detection_threshold )
                {
                    SignalBehavior.setSingleBehavior( pCell, "migration speed", 0 );
                }
                SignalBehavior.setSingleBehavior( pCell, "migration bias", unattached_worker_migration_bias );
                SignalBehavior.setSingleBehavior( pCell, "chemotactic response to oxygen", 0 );
                SignalBehavior.setSingleBehavior( pCell, "chemotactic response to chemoattractant", 1 );
            }
        }
    }

    public static class TherapyEvent extends Event
    {
        public TherapyEvent(double executionTime)
        {
            super( executionTime );
        }

        @Override
        public void execute(Model model) throws Exception
        {
            System.out.println( "Therapy started!" );
            model.setSaveInterval( model.getParameterDouble( "save_interval_after_therapy_start" ) ); // 3.0; 
            introduce_biorobots( model );
        }
    }
}