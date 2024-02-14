package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    private static Map<Integer, Integer> typeToIndex = new HashMap<>();
    private static Map<String, CellDefinition> cellDefinitionNames = new HashMap<>();
    private static Map<Integer, CellDefinition> cellDefinitionTypes = new HashMap<>();

    public int type;
    public String name;

    public boolean isMovable;

    private Microenvironment m;
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
        cellDefinitionNames.put( cd.name, cd );
        cellDefinitionTypes.put( cd.type, cd );
        typeToIndex.put( cd.type, cellDefinitions.size() );
        cellDefinitions.add( cd );
        sync();
    }

    public static void clearCellDefinitions()
    {
        cellDefinitions.clear();
        cellDefinitionNames.clear();
        cellDefinitionTypes.clear();
        typeToIndex.clear();
        sync();
    }

    public static Iterable<CellDefinition> getCellDefinitions()
    {
        return cellDefinitions;
    }

    public static CellDefinition getCellDefinitionByIndex(int index)
    {
        return cellDefinitions.get( index );
    }

    public static int getCellDefinitionIndex(int type)
    {
        return typeToIndex.get( type );
    }

    public static CellDefinition getCellDefinition(int type)
    {
        return cellDefinitionTypes.get( type );
    }

    public static Set<String> getCellDefinitionNames()
    {
        return cellDefinitionNames.keySet();
    }

    public static CellDefinition getCellDefinition(String name)
    {
        return cellDefinitionNames.get( name );
    }

    private static void sync()
    {
        for( CellDefinition cd : cellDefinitions )
        {
            cd.phenotype.sync();
        }
    }

    public CellDefinition()
    {
        isMovable = true;
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
        functions.set_orientation = null;
    }

    public CellDefinition(Microenvironment m, int type, String name)
    {
        this();
        this.m = m;
        this.type = type;
        this.name = name;
        phenotype.sync( m );
    }

    public static int findCellDefinitionIndex(String name)
    {
        CellDefinition cd = cellDefinitionNames.get( name );
        if( cd != null )
            return cd.type;
        return -1;
    }

    public CellDefinition clone(String name, int type, Microenvironment m)
    {
        CellDefinition result = new CellDefinition( m, type, name );
        result.isMovable = this.isMovable;
        result.parameters = parameters.clone();
        result.custom_data = custom_data.clone();
        result.phenotype = phenotype.clone();
        result.functions = functions.clone();
        result.phenotype.sync( m );
        return result;
    }

    @Override
    public String toString()
    {
        return name + " ( " + type + " ) ";
    }

    public Microenvironment getMicroenvironment()
    {
        return m;
    }

    //    String display(  )
    //    {
    //        StringBuilder sb = new StringBuilder();
    //            sb.append( " type:" + type + " name: " + name ); 
    //
    //            // summarize cycle model 
    //            if( phenotype.cycle!= null )
    //            {
    //                sb.append( "\t cycle model: " + phenotype.cycle.name  
    //                    + " (code=" + phenotype.cycle.code + ")" ); 
    //                    
    //                // let's show the transition rates 
    //                CycleModel pCM = (phenotype.cycle ); 
    //                CycleData pCMD = (phenotype.cycle.data ); 
    //                for( int n=0 ; n < pCM.phases.size() ; n++ )
    //                {
    //                    sb.append( "\t\tPhase " + n + ": " + pCM.phases.get(n).name+"\n" ); 
    //                }
    //                sb.append( "\t\tCycle transitions: \n" 
    //                   + "\t\t-----------------\n" ); 
    //                for( int n=0 ; n < pCM.phase_links.size() ; n++ )
    //                {
    //                    for( int k=0; k < pCM.phase_links.get(n).size() ; k++ )
    //                    {
    //                        int start = pCM.phase_links[n][k].start_phase_index;
    //                        int end = pCM.phase_links[n][k].end_phase_index; 
    //                        sb.append( "\t\t" + pCM.phases[start].name + " -. " 
    //                            + pCM.phases[end].name + " w mean duration " 
    //                            + 1.0 / pCMD.transition_rate( start,end) + " min" ); 
    //                    }
    //                }           
    //                
    //            }
    //            else
    //            {   sb.append( "\t cycle model not initialized" ); } 
    //
    //            // summarize death models 
    //            sb.append( "\t death models: " ); 
    //            for( int k=0 ; k < phenotype.death.models.size(); k++ )
    //            {
    //                sb.append( "\t\t" + k + " : " + phenotype.death.models[k].name 
    //                + " (code=" + phenotype.death.models[k].code + ")" 
    //                + " with rate " + phenotype.death.rates[k] + " 1/min" ); 
    //
    //                CycleModel pCM = (phenotype.death.models[k] ); 
    //                CycleData pCMD = (phenotype.death.models[k].data ); 
    //
    //                
    //                sb.append( "\t\tdeath phase transitions: \n" ) 
    //                   + "\t\t------------------------" ); 
    //                for( int n=0 ; n < pCM.phase_links.size() ; n++ )
    //                {
    //                    for( int k=0; k < pCM.phase_links[n].size() ; k++ )
    //                    {
    //                        int start = pCM.phase_links[n][k].start_phase_index;
    //                        int end = pCM.phase_links[n][k].end_phase_index; 
    //                        sb.append( "\t\t" + pCM.phases[start].name + " -. " 
    //                            + pCM.phases[end].name + " w mean duration " 
    //                            + 1.0 / pCMD.transition_rate( start,end) + " min" ); 
    //                    }
    //                }           
    //                
    //                
    //                
    //            }
    //            
    //            // summarize functions 
    //            CellFunctions pCF = (functions); 
    //            sb.append( "\t key functions: " ); 
    //            sb.append( "\t\t migration bias rule: "; display_ptr_as_bool( pCF.update_migration_bias  ); 
    //            sb.append( "\n"); 
    //            sb.append( "\t\t custom rule: "; display_ptr_as_bool( pCF.custom_cell_rule  ); 
    //            sb.append( "\n");
    //            sb.append( "\t\t phenotype rule: "; display_ptr_as_bool( pCF.update_phenotype ); 
    //            sb.append( "\n"); 
    //            sb.append( "\t\t volume update function: "; display_ptr_as_bool( pCF.volume_update_function  ); 
    //            sb.append( "\n"); 
    //            sb.append( "\t\t mechanics function: "; display_ptr_as_bool( pCF.update_velocity  ); 
    //            sb.append( "\n");
    //            sb.append( "\t\t contact function: "; display_ptr_as_bool( pCF.contact_function ); 
    //            sb.append( "\n"); 
    //            
    //            // summarize motility 
    //            
    //            Motility pM = (phenotype.motility); 
    //            String val = "true";
    //            if( pM.isMotile == false )
    //            { val = "false"; } 
    //        
    //            String dimen = "3D"; 
    //            if( pM.restrictTo2D == true )
    //            { dimen = "2D"; } 
    //
    //            sb.append( "\tmotility (enabled: " + val + " in " + dimen + ")\n" 
    //                + "\t\tspeed: " + pM.migrationSpeed + " micron/min\n"  
    //                + "\t\tbias: " + pM.migrationBias + " \n" 
    //                + "\t\tpersistence time: " + pM.persistenceTime + " min\n"
    //                + "\t\tchemotaxis (enabled: ");
    //                
    //                val = "true" ;
    //                if( functions.update_migration_bias != chemotaxis_function )
    //                { val = "false"; } 
    //            sb.append( val + ")\n"); 
    //                + "\t\t\talong " 
    //                + pM.chemotaxis_direction + " * grad(" 
    //                + microenvironment.density_names[ pM.chemotaxis_index ] + ") " ); 
    //                
    //            // secretion
    //            
    //            
    //            
    //            // mechanics
    //
    //            Mechanics pMech = (phenotype.mechanics); 
    //
    //            sb.append( "\tmechanics:\n" 
    //                + "\t\tcell_cell_adhesion_strength: " + pMech.cellCellAdhesionStrength + "\n" 
    //                + "\t\tcell_cell_repulsion_strength: " + pMech.cellCellRepulsionStrength + "\n" 
    //                + "\t\trel max adhesion dist: " + pMech.relMaxAdhesionDistance + "\n" 
    //                + "\t\tcell_BM_adhesion_strength: " + pMech.cellBMAdhesionStrength + "\n" 
    //                + "\t\tcell_BM_repulsion_strength: " + pMech.cellBMRepulsionStrength + "\n" 
    //                + "\t\tattachment_elastic_constant: " + pMech.attachmentElasticConstant + "\n" 
    //                + "\t\tattachment_rate: " + pMech.attachmentRate + "\n" 
    //                + "\t\tdetachment_rate: " + pMech.detachmentRate );
    //            
    //            // size 
    //        
    //
    //            // intracellular
    //            if (phenotype.intracellular != null)
    //            {
    //                phenotype.intracellular.display();
    //            }
    //            
    //            CustomCellData pCCD = (custom_data); 
    //            sb.append( "\tcustom data: " ); 
    //            for( int k=0; k < pCCD.variables.size(); k++)
    //            {
    //                sb.append( "\t\t" + pCCD.variables.get(k) ); 
    //            }
    //            sb.append( "\tcustom vector data: " ); 
    //            for( int k=0; k < pCCD.vectorVariables.size(); k++)
    //            {
    //                sb.append( "\t\t" + pCCD.vectorVariables.get(k) ); 
    //            }
    //            sb.append( "\t\t\tNOTE: custom vector data will eventually be merged with custom data" );
    //            return sb.toString();
    //    }
}