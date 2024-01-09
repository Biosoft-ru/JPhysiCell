package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ru.biosoft.physicell.biofvm.Microenvironment;

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
public class CellDefinition
{
    //Static, move to registry of some kind
    private static List<CellDefinition> cellDefinitions = new ArrayList<>();
    private static Map<String, CellDefinition> cellDefinitionsMap = new HashMap<>();

    public int type;
    public String name;

    boolean is_movable;

    public CellParameters parameters = new CellParameters();
    public CustomCellData custom_data = new CustomCellData();
    public CellFunctions functions = new CellFunctions();
    public Phenotype phenotype = new Phenotype();


    public static int getDefinitionsCount()
    {
        return cellDefinitions.size();
    }

    public static void registerCellDefinition(CellDefinition cd)
    {
        cd.type = cellDefinitions.size();
        cellDefinitions.add( cd );
        cellDefinitionsMap.put( cd.name, cd );
        sync();
    }

    public static void clearCellDefinitions()
    {
        cellDefinitions.clear();
        cellDefinitionsMap.clear();
        sync();
    }

    public static CellDefinition getCellDefinition(int type)
    {
        return cellDefinitions.get( type );
    }

    public static CellDefinition getCellDefinition(String name)
    {
        return cellDefinitionsMap.get( name );
    }

    private static void sync()
    {
        for( CellDefinition cd : cellDefinitions )
        {
            cd.phenotype.sync();
        }
    }

    public CellDefinition(Microenvironment m, String name)
    {
        // set up the default parameters, the default Cell_Parameters constructor should take care of this
        this.type = -1; //not registered
        this.name = name;
        is_movable = true;
        //        parameters.pReference_live_phenotype = phenotype; //TODO: check

        // set up the default custom data, the default Custom_Cell_Data constructor should take care of this

        // set up the default functions 
        functions.instantiate_cell = null;
        functions.updateVolume = null; // standard_volume_update_function;
        functions.update_migration_bias = null;

        functions.updatePhenotype = null;
        functions.custom_cell_rule = null;

        functions.updateVelocity = null; // standard_update_cell_velocity;
        functions.add_cell_basement_membrane_interactions = null;
        functions.calculate_distance_to_membrane = null;
        //
        //        // bug fix July 31 2023
        //                functions.plot_agent_SVG = standard_agent_SVG;
        //                functions.plot_agent_legend = standard_agent_legend;
        //        // bug fix July 31 2023
        //        
        functions.set_orientation = null;

        phenotype.sync( m );
        // new March 2022 : make sure Cell_Interactions and Cell_Transformations are appropriately sized. Same on motiltiy. 
        //        phenotype.sync();
    }

    //    public CellDefinition( CellDefinition cd )
    //    {
    //        // set the microenvironment pointer 
    //        pMicroenvironment = cd.pMicroenvironment;
    //
    //        // set up the default parameters 
    //            // the default Cell_Parameters constructor should take care of this
    //            
    //        type = cd.type; 
    //        name = cd.name; 
    //         
    //        parameters = cd.parameters;
    //        custom_data = cd.custom_data; 
    //        functions = cd.functions; 
    //        phenotype = cd.phenotype; 
    //        
    //        // this is the whole reason we need ot make a copy constructor 
    //        parameters.pReference_live_phenotype = &phenotype; 
    //        
    //        cell_definitions_by_index.push_back( this ); 
    //        
    //        return; 
    //    }
    //
    //    CellDefinition::operator=(const Cell_Definition&cd)
    //    {
    //        // set the microenvironment pointer 
    //        pMicroenvironment = cd.pMicroenvironment;
    //
    //        // set up the default parameters 
    //            // the default Cell_Parameters constructor should take care of this
    //            
    //        type = cd.type; 
    //        name = cd.name; 
    //         
    //        parameters = cd.parameters;
    //        custom_data = cd.custom_data; 
    //        functions = cd.functions; 
    //        phenotype = cd.phenotype; 
    //        
    //        // this is the whole reason we need ot make a copy constructor 
    //        parameters.pReference_live_phenotype = &phenotype; 
    //        
    //        // commented out on March 10, 2020 
    //        // cell_definitions_by_index.push_back( this ); 
    //        
    //        return *this; 
    //    }
}
