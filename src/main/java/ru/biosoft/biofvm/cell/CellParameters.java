package ru.biosoft.biofvm.cell;

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
public class CellParameters implements Cloneable
{
    // oxygen values (in mmHg) for critical phenotype changes
    double o2_hypoxic_threshold; // value at which hypoxic signaling starts
    double o2_hypoxic_response; // value at which omics changes are observed 
    double o2_hypoxic_saturation; // value at which hypoxic signalign saturates 
    // o2_hypoxic_saturation < o2_hypoxic_threshold

    double o2_proliferation_saturation; // value at which extra o2 does not increase proliferation
    double o2_proliferation_threshold; // value at which o2 is sufficient for proliferation

    double o2_reference; // physioxic reference value, in the linked reference Phenotype
    // o2_proliferation_threshold < o2_reference < o2_proliferation_saturation; 

    double o2_necrosis_threshold; // value at which cells start experiencing necrotic death 
    double o2_necrosis_max; // value at which necrosis reaches its maximum rate 
    // o2_necrosis_max < o2_necrosis_threshold

    Phenotype pReference_live_phenotype; // reference live phenotype (typically physioxic) 
    Phenotype pReference_necrotic_phenotype; // reference live phenotype (typically physioxic) 

    // necrosis parameters (may evenually be moved into a reference necrotic phenotype 
    public double max_necrosis_rate; // deprecate
    int necrosis_type; // deprecate 

    public CellParameters()
    {
        o2_hypoxic_threshold = 15.0; // HIF-1alpha at half-max around 1.5-2%, and tumors often are below 2%
        o2_hypoxic_response = 8.0; // genomic / proteomic changes observed at 7-8 mmHg 
        o2_hypoxic_saturation = 4.0; // maximum HIF-1alpha at 0.5% o2 (McKeown)

        o2_necrosis_threshold = 5.0;
        o2_necrosis_max = 2.5;

        o2_proliferation_threshold = 5.0; // assume no proliferation at same level as starting necrosis 
        o2_proliferation_saturation = 160.0; // 5% = 38, 21% = 160 mmHg 
        o2_reference = 160.0; // assume all was measured in normoxic 21% o2 

        //        pReference_live_phenotype = NULL; // reference live (usually physioxic) phenotype //TODO: uncomment 

        // necrosis parameters 

        max_necrosis_rate = 1.0 / ( 6.0 * 60.0 ); // assume cells survive 6 hours in very low oxygen 
        //        necrosis_type = PhysiCell_constants::deterministic_necrosis;;//TODO: uncomment

        return;
    }

    @Override
    public CellParameters clone()
    {
        try
        {
            return (CellParameters)super.clone();
        }
        catch( CloneNotSupportedException e )
        {
            e.printStackTrace();
            return null;
        }
    }
}
