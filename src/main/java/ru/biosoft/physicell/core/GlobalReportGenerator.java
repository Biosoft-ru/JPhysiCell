package ru.biosoft.physicell.core;

public class GlobalReportGenerator
{
    public String getReportElements(Model model) throws Exception
    {
        return PhysiCellUtilities.getCurrentTime() + "\tElapsed\t" + ( System.currentTimeMillis() - model.startTime ) / 1000 + "\tTime:\t"
                + (int)Math.round( model.curTime ) + "\tCells\t" + model.getMicroenvironment().getAgentsCount();
    }
}
