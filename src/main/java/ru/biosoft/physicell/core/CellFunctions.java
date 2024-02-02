package ru.biosoft.physicell.core;

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
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
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
public class CellFunctions
{
    //    public CycleModel cycleModel; TODO: temporarily commented

    public instantiate_cell instantiate_cell;
    public volume_update_function updateVolume;
    public update_migration_bias update_migration_bias;
    public custom_cell_rule custom_cell_rule;
    public update_phenotype updatePhenotype;
    public pre_update_intracellular pre_update_intracellular;
    public post_update_intracellular post_update_intracellular;
    public update_velocity updateVelocity;
    public add_cell_basement_membrane_interactions add_cell_basement_membrane_interactions;
    public calculate_distance_to_membrane calculate_distance_to_membrane;
    public set_orientation set_orientation;
    public contact_function contact_function;

    /* prototyping / beta in 1.5.0 */
    /*  
        void (*internal_substrate_function)(Cell* pCell, Phenotype& phenotype , double dt ); 
        void (*molecular_model_function)(Cell* pCell, Phenotype& phenotype , double dt ); 
    */


    //           void (*plot_agent_SVG)(std::ofstream& os, Cell* pCell, double z_slice, std::vector<std::string> (*cell_coloring_function)(Cell*), double X_lower, double Y_lower);
    //           void (*plot_agent_legend)(std::ofstream& os, Cell_Definition* cell_def, double& cursor_x, double& cursor_y, std::vector<std::string> (*cell_coloring_function)(Cell*), double temp_cell_radius);

    public CellFunctions clone()
    {

        CellFunctions result = new CellFunctions();
        try
        {
            result.instantiate_cell = instantiate_cell == null ? null : instantiate_cell.getClass().newInstance();
            result.updateVolume = updateVolume == null ? null : updateVolume.getClass().newInstance();
            result.update_migration_bias = update_migration_bias == null ? null : update_migration_bias.getClass().newInstance();
            result.updatePhenotype = updatePhenotype == null ? null : updatePhenotype.getClass().newInstance();
            result.pre_update_intracellular = pre_update_intracellular == null ? null : pre_update_intracellular.getClass().newInstance();
            result.post_update_intracellular = post_update_intracellular == null ? null
                    : post_update_intracellular.getClass().newInstance();
            result.updateVelocity = updateVelocity == null ? null : updateVelocity.getClass().newInstance();
            result.add_cell_basement_membrane_interactions = add_cell_basement_membrane_interactions == null ? null
                    : add_cell_basement_membrane_interactions.getClass().newInstance();
            result.calculate_distance_to_membrane = calculate_distance_to_membrane == null ? null
                    : calculate_distance_to_membrane.getClass().newInstance();
            result.set_orientation = set_orientation == null ? null : set_orientation.getClass().newInstance();
            result.contact_function = contact_function == null ? null : contact_function.getClass().newInstance();
            result.custom_cell_rule = custom_cell_rule == null ? null : custom_cell_rule.getClass().newInstance();
        }
        catch( Exception ex )
        {
            ex.printStackTrace();
        }
        return result;
    }

    @FunctionalInterface
    public static interface volume_update_function
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface update_migration_bias
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface update_phenotype
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception;
    }

    @FunctionalInterface
    public static interface update_velocity
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface custom_cell_rule
    {
		public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception;
    }

    @FunctionalInterface
    public static interface pre_update_intracellular
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface post_update_intracellular
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface add_cell_basement_membrane_interactions
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface calculate_distance_to_membrane
    {
        public double execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface set_orientation
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt);
    }

    @FunctionalInterface
    public static interface contact_function
    {
		public void execute(Cell pCell, Phenotype phenotype, Cell cell2, Phenotype phenotype2, double dt);
    }

    @FunctionalInterface
    public static interface instantiate_cell
    {
        public Cell execute();
    }
}