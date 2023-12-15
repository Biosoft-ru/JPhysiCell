package ru.biosoft.biofvm.cell;

import ru.biosoft.biofvm.VectorUtil;

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
public class CellInteractions implements Cloneable
{
    // phagocytosis parameters (e.g., macrophages)
    double dead_phagocytosis_rate;

    double[] live_phagocytosis_rates;
    // attack parameters (e.g., T cells)

    double[] attack_rates;
    // do I attack cell type j? 

    double[] immunogenicities; // new! 
    // how immnogenic am I to cell type j? 

    double damage_rate;
    // cell fusion parameters 
    double[] fusion_rates;

    public CellInteractions()
    {
        dead_phagocytosis_rate = 0.0;
        live_phagocytosis_rates = new double[] {0.0};
        damage_rate = 1.0;
        attack_rates = new double[] {0.0};
        immunogenicities = new double[] {1};
        fusion_rates = new double[] {0.0};
    }

    void sync_to_cell_definitions()
    {
        //        extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
        int number_of_cell_defs = Cell.cell_definition_indices_by_name.size();
        
        if( live_phagocytosis_rates.length != number_of_cell_defs )
        {
            live_phagocytosis_rates = VectorUtil.resize( live_phagocytosis_rates, number_of_cell_defs );
            attack_rates = VectorUtil.resize( attack_rates, number_of_cell_defs );
            fusion_rates = VectorUtil.resize( fusion_rates, number_of_cell_defs );
            immunogenicities = VectorUtil.resize( immunogenicities, number_of_cell_defs );

            //            live_phagocytosis_rates.resize( number_of_cell_defs, 0.0); 
            //            attack_rates.resize( number_of_cell_defs, 0.0); 
            //            fusion_rates.resize( number_of_cell_defs, 0.0); 
            //            immunogenicities.resize( number_of_cell_defs , 1.0 ); 
        }
    }
    public double live_phagocytosis_rate(String type_name)
    {
        int n = Cell.cell_definition_indices_by_name.get( type_name );
        return live_phagocytosis_rates[n];
    }
    public double attack_rate(String type_name)
    {
        int n = Cell.cell_definition_indices_by_name.get( type_name );
        return attack_rates[n];
    }

    public double fusion_rate(String type_name)
    {
        int n = Cell.cell_definition_indices_by_name.get( type_name );
        return fusion_rates[n];
    }

    public double immunogenicity(String type_name)
    {
        int n = Cell.cell_definition_indices_by_name.get( type_name );
        return immunogenicities[n];
    }

    // ease of access 
    //    double&Cell_Interactions::live_phagocytosis_rate( std::string type_name )
    //    {
    //        extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
    //        int n = cell_definition_indices_by_name[type_name]; 
    //        // std::cout << type_name << " " << n << std::endl; 
    //        return live_phagocytosis_rates[n]; 
    //    }
    //
    //    double& Cell_Interactions::attack_rate( std::string type_name ) 
    //    {
    //        extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
    //        int n = cell_definition_indices_by_name[type_name]; 
    //        return attack_rates[n]; 
    //    }
    //
    //    double& Cell_Interactions::fusion_rate( std::string type_name )
    //    {
    //        extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
    //        int n = cell_definition_indices_by_name[type_name]; 
    //        return fusion_rates[n]; 
    //    }
    //
    //    double& Cell_Interactions::immunogenicity( std::string type_name )
    //    {
    //        extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
    //        int n = cell_definition_indices_by_name[type_name]; 
    //        return immunogenicities[n]; 
    //    }
    @Override
    public CellInteractions clone()
    {
        try
        {
            CellInteractions result = (CellInteractions)super.clone();
            result.live_phagocytosis_rates = this.live_phagocytosis_rates.clone();
            result.attack_rates = this.attack_rates.clone();
            result.immunogenicities = this.immunogenicities.clone();
            result.fusion_rates = this.fusion_rates.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}