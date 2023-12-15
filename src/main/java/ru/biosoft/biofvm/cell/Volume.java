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
public class Volume implements Cloneable
{
    double total;
    double solid;
    double fluid;
    double fluid_fraction;

    double nuclear;
    double nuclear_fluid;
    double nuclear_solid;

    double cytoplasmic;
    double cytoplasmic_fluid;
    double cytoplasmic_solid;

    double calcified_fraction;

    double cytoplasmic_to_nuclear_ratio;

    double rupture_volume; // in volume units 

    //
    // parameters that can be set by users 
    //
    double cytoplasmic_biomass_change_rate;
    double nuclear_biomass_change_rate;
    double fluid_change_rate;

    double calcification_rate;

    double target_solid_cytoplasmic;
    double target_solid_nuclear;
    double target_fluid_fraction;

    double target_cytoplasmic_to_nuclear_ratio;

    public double[] test = new double[] {3, 4};
    public double relative_rupture_volume;
    // the volume ratio (compared to initial volume at time of death) 
    // at which a cell ruptures / lyses / bursts. 

    public Volume()
    {
        // reference parameter values for MCF-7, in cubic microns 
        fluid_fraction = 0.75;

        total = 2494;
        fluid = fluid_fraction * total;
        solid = total - fluid;

        nuclear = 540.0;

        nuclear_fluid = fluid_fraction * nuclear;
        nuclear_solid = nuclear - nuclear_fluid;

        cytoplasmic = total - nuclear;
        cytoplasmic_fluid = fluid_fraction * cytoplasmic;
        cytoplasmic_solid = cytoplasmic - cytoplasmic_fluid;

        // rates are in units of 1/min 
        cytoplasmic_biomass_change_rate = 0.27 / 60.0;
        nuclear_biomass_change_rate = 0.33 / 60.0;
        fluid_change_rate = 3.0 / 60.0;

        calcified_fraction = 0.0;

        calcification_rate = 0.0;

        target_solid_cytoplasmic = cytoplasmic_solid;
        target_solid_nuclear = nuclear_solid;
        target_fluid_fraction = fluid_fraction;

        cytoplasmic_to_nuclear_ratio = cytoplasmic / ( 1e-16 + nuclear );
        target_cytoplasmic_to_nuclear_ratio = cytoplasmic_to_nuclear_ratio;

        // the cell bursts at these volumes 
        relative_rupture_volume = 2.0;
        // as fraction of volume at entry to the current phase
        rupture_volume = relative_rupture_volume * total; // in volume units 

    }

    public void multiply_by_ratio(double ratio)
    {
        total *= ratio;
        solid *= ratio;
        fluid *= ratio;

        nuclear *= ratio;
        nuclear_fluid *= ratio;
        nuclear_solid *= ratio;

        cytoplasmic *= ratio;
        cytoplasmic_fluid *= ratio;
        cytoplasmic_solid *= ratio;

        rupture_volume *= ratio;

        target_solid_nuclear *= ratio;
        target_solid_cytoplasmic *= ratio;
    }

    public void divide()
    {
        multiply_by_ratio( 0.5 );
    }

    @Override
    public Volume clone()
    {
        try
        {
            return (Volume)super.clone();
        }
        catch( CloneNotSupportedException e )
        {
            throw new InternalError( e );
        }
    }
}