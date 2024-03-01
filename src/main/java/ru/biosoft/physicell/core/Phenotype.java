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
    boolean flaggedForDivision;
    boolean flaggedForRemoval;

    public CycleModel cycle = new CycleModel();
    public Death death = new Death();
    public Volume volume = new Volume();
    public Geometry geometry = new Geometry();
    public Mechanics mechanics = new Mechanics();
    public Motility motility = new Motility();
    public Secretion secretion = new Secretion();
    public Molecular molecular = new Molecular();
    public CellInteractions cellInteractions = new CellInteractions();
    public CellTransformations cellTransformations = new CellTransformations();
    public Intracellular intracellular = null;//new Intracellular();

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
    }

    public void initialize(int cellDefinitionSize)
    {
        cellInteractions.initialize( cellDefinitionSize );
        cellTransformations.initialize( cellDefinitionSize );
        mechanics.initialize( cellDefinitionSize );
    }

    public Phenotype()
    {
        flaggedForDivision = false;
        flaggedForRemoval = false;
        intracellular = null; // sync the molecular stuff here automatically? 
    }

    @Override
    public Phenotype clone()
    {
        Phenotype result = new Phenotype();
        result.flaggedForDivision = flaggedForDivision;
        result.flaggedForRemoval = flaggedForRemoval;
        result.cycle = cycle.clone();
        result.death = death.clone();
        result.volume = volume.clone();
        result.geometry = geometry.clone();
        result.mechanics = mechanics.clone();
        result.motility = motility.clone();
        result.secretion = secretion.clone();
        result.molecular = molecular.clone();
        result.cellInteractions = cellInteractions.clone();
        result.cellTransformations = cellTransformations.clone();
        if( intracellular != null )
            result.intracellular = intracellular.clone();
        return result;
    }
}