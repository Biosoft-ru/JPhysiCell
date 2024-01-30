package ru.biosoft.physicell.core;

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
public class Phenotype implements Cloneable
{
    boolean flagged_for_division;
    boolean flagged_for_removal;

    public CycleModel cycle = new CycleModel();
    public Death death = new Death();
    public Volume volume = new Volume();
    public Geometry geometry = new Geometry();
    public Mechanics mechanics = new Mechanics();
    public Motility motility = new Motility();
    public Secretion secretion = new Secretion();
    public Molecular molecular = new Molecular();
    public CellInteractions cell_interactions = new CellInteractions();
    public CellTransformations cell_transformations = new CellTransformations();

    // We need it to be a pointer to allow polymorphism
    // then this object could be a MaBoSSIntracellular, or a RoadRunnerIntracellular
    public Intracellular intracellular = new Intracellular();

    //    public void sync(CellFunctions functions)
    //    {
    //        cycle = functions.cycleModel;
    //    }

    /**
     * Synchronize all parts with microenvironment (densities)  
     */
    public void sync(Microenvironment m)
    {
        motility.sync( m );
        secretion.sync( m );
        molecular.sync( m );
    }

    /**
     * Synchronize all parts with Cell Definition
     */
    public void sync()
    {
        this.initialize( CellDefinition.getDefinitionsCount() );
        //        cell_interactions.sync_to_cell_definitions();
        //        cell_transformations.sync_to_cell_definitions();
        //        mechanics.sync_to_cell_definitions();
    }

    public void initialize(int cellDefinitionSize)
    {
        cell_interactions.initialize( cellDefinitionSize );
        cell_transformations.initialize( cellDefinitionSize );
        mechanics.initialize( cellDefinitionSize );
    }

    public Phenotype()
    {
        flagged_for_division = false;
        flagged_for_removal = false;

        // sync the molecular stuff here automatically? 
        intracellular = null;
    }

    @Override
    public Phenotype clone()
    {
        Phenotype result = new Phenotype();
        result.flagged_for_division = flagged_for_division;
        result.flagged_for_removal = flagged_for_removal;
        result.cycle = cycle.clone();
        result.death = death.clone();
        result.volume = volume.clone();
        result.geometry = geometry.clone();
        result.mechanics = mechanics.clone();
        result.motility = motility.clone();
        result.secretion = secretion.clone();
        result.molecular = molecular.clone();
        result.cell_interactions = cell_interactions.clone();
        result.cell_transformations = cell_transformations.clone();
        if( intracellular != null )
        result.intracellular = intracellular.clone();//
        //        flagged_for_division = p.flagged_for_division;
        //        flagged_for_removal = p.flagged_for_removal;
        //        cycle = p.cycle;
        //        death = p.death;
        //        volume = p.volume;
        //        geometry = p.geometry;
        //        mechanics = p.mechanics;
        //        motility = p.motility;
        //        secretion = p.secretion;

        //        molecular = p.molecular;
        //        
        //        delete intracellular;
        //        
        //        if (p.intracellular != NULL)
        //        { intracellular = p.intracellular->clone(); }
        //        else
        //        { intracellular = NULL; }
        //        
        //        cell_interactions = p.cell_interactions; 
        //        cell_transformations = p.cell_transformations; 
        return result;
    }
}
