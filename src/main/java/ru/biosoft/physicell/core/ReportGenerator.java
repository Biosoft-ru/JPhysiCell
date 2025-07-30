package ru.biosoft.physicell.core;

public class ReportGenerator
{
    public String[] getReportHeaderElements()
    {
        return new String[] {"X", "Y", "Z", "Cycle", "Elapsed"};
    }

    public Object[] getReportElements(Cell cell) throws Exception
    {
        return new Object[] {cell.position[0], cell.position[1], cell.position[2], cell.phenotype.cycle.currentPhase().name,
                cell.phenotype.cycle.data.elapsedTimePhase};
    }
}
