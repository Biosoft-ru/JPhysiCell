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
public class Intracellular implements Cloneable
{

    String intracellular_type; // specified in XML <intracellular type="...">:  "maboss", "sbml", ...
    // bool enabled; 

    // ==========  specific to SBML ==============
    // std::string sbml_filename;


    // ================  generic  ================
    // This function parse the xml cell definition
    //        virtual void initialize_intracellular_from_pugixml(pugi::xml_node& node) = 0;
    //        
    //        // This function initialize the model, needs to be called on each cell once created
    //        virtual void start() = 0;
    //        
    //        // This function checks if it's time to update the model
    //        virtual bool need_update() = 0;

    public void start()
    {

    }

    public boolean need_update()
    {
        return false;
    }

    public void update()
    {

    }

    public void update(Cell cell, Phenotype phenotype, double dt)
    {

    }

    public void inherit(Cell cell)
    {

    }
    //
    //        // This function update the model for the time_step defined in the xml definition
    //        virtual void update() = 0;
    //        virtual void update(Cell* cell, Phenotype& phenotype, double dt) = 0;
    //
    //        // This function deals with inheritance from mother to daughter cells
    //        virtual void inherit(Cell* cell) = 0;
    //
    //        // Get value for model parameter
    //        virtual double get_parameter_value(std::string name) = 0;
    //        
    //        // Set value for model parameter
    //        virtual void set_parameter_value(std::string name, double value) = 0;
    //
    //        virtual std::string get_state() = 0;  
    //        
    //        virtual void display(std::ostream& os) = 0;
    //
    //        virtual Intracellular* clone() = 0;
    //        
    //        virtual ~Intracellular(){};
    //        
    //
    //        // ================  specific to "maboss" ================
    //        virtual bool has_variable(std::string name) = 0; 
    //        virtual bool get_boolean_variable_value(std::string name) = 0;
    //        virtual void set_boolean_variable_value(std::string name, bool value) = 0;
    //        // virtual bool get_double_variable_value(std::string name) = 0;
    //        // virtual void set_double_variable_value(std::string name, bool value) = 0;
    //        virtual void print_current_nodes() = 0;
    //        
    //
    //        // ================  specific to "roadrunner" ================
    //        virtual int update_phenotype_parameters(PhysiCell::Phenotype& phenotype) = 0;
    //        virtual int validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype) = 0;
    //        virtual int validate_SBML_species() = 0;
    //        virtual int create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype) = 0;
    @Override
    public Intracellular clone()
    {
        try
        {
            return (Intracellular)super.clone();
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}
