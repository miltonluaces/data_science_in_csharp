﻿#region Imports

using System;
using System.Data;
using System.Drawing;
using System.Collections;
using System.Text;
using System.IO;

#endregion

namespace AibuSet {

    internal class Reporting {
    
        #region Fields
        
        private DataSet reportSource;
        private ArrayList sections;
        private string reportTitle;
        private string newline;
        private ArrayList reportFields;
        private int iLevel = 0;
        private string gradientStyle;
        private StringBuilder htmlContent;
        internal string ReportFont;
        private Hashtable totalList;
        internal ArrayList TotalFields;
        internal bool IncludeTotal;
        //Chart fields
        internal bool IncludeChart;
        internal string ChartTitle;
        internal bool ChartShowAtBottom;
        internal string ChartChangeOnField;
        internal string ChartValueField = "Count";
        internal bool ChartShowBorder;
        internal string ChartLabelHeader = "Label";
        internal string ChartPercentageHeader = "Percentage";
        internal string ChartValueHeader = "Value";
        
        #endregion

        #region Constructor

        internal Reporting() {
            htmlContent = new StringBuilder();
            newline = "\n";
            sections = new ArrayList();
            reportFields = new ArrayList();
            ReportFont = "Arial";
            gradientStyle = "FILTER: progid:DXImageTransform.Microsoft.Gradient(gradientType=1,startColorStr=BackColor,endColorStr=#ffffff)";
            totalList = new Hashtable();
            TotalFields = new ArrayList();
        }

        #endregion

        #region Properties

        internal DataSet ReportSource {
            get { return reportSource; }
            set { reportSource = value; }
        }

        internal ArrayList Sections {
            get { return sections; }
            set { sections = value; }
        }

        internal string ReportTitle {
            get { return reportTitle; }
            set { reportTitle = value; }
        }

        internal ArrayList ReportFields {
            get { return reportFields; }
            set { reportFields = value; }
        }
        
        #endregion

        #region internal Methods

        internal string GenerateReport() {
            foreach (Field fld in this.ReportFields) {
                if (!this.TotalFields.Contains(fld.FieldName) && fld.isTotalField) { TotalFields.Add(fld.FieldName); }
            }
            WriteTitle();
            WriteSections();
            WriteFooter();
            return htmlContent.ToString();
        }

        internal bool Save(string fileName) {
            try {
                GenerateReport();
                StreamWriter sw = new StreamWriter(fileName);
                sw.Write(htmlContent.ToString());
                sw.Flush();
                sw.Close();
                return true;
            }
            catch (Exception exp) {
                System.Diagnostics.Debug.WriteLine(exp.Message);
                return false;
            }
        }

        #endregion

        #region Private Methods
        
        private void WriteTitle() {
            htmlContent.Append("<HTML><HEAD><TITLE>Report - " + reportTitle + "</TITLE></HEAD>" + newline);
            htmlContent.Append("<STYLE>" + newline);
            htmlContent.Append(" .TableStyle { border-collapse: collapse } " + newline);
            htmlContent.Append(" .TitleStyle { font-family: " + ReportFont + "; font-size:15pt } " + newline);
            htmlContent.Append(" .SectionHeader {font-family: " + ReportFont + "; font-size:10pt } " + newline);
            htmlContent.Append(" .DetailHeader {font-family: " + ReportFont + "; font-size:9pt } " + newline);
            htmlContent.Append(" .DetailData  {font-family: " + ReportFont + "; font-size:9pt } " + newline);
            htmlContent.Append(" .ColumnHeaderStyle  {font-family: " + ReportFont + "; font-size:9pt; border-style:outset; border-width:1} " + newline);
            htmlContent.Append("</STYLE>" + newline);
            htmlContent.Append("<BODY TOPMARGIN=0 LEFTMARGIN=0 RIGHTMARGIN=0 BOTTOMMARGIN=0>" + newline);
            htmlContent.Append("<TABLE Width='100%' style='FILTER: progid:DXImageTransform.Microsoft.Gradient(gradientType=1,startColorStr=#a9d4ff,endColorStr=#ffffff)' Cellpadding=5><TR><TD>");
            htmlContent.Append("<font face='" + ReportFont + "' size=6>" + ReportTitle + "</font>");
            htmlContent.Append("</TD></TR></TABLE>" + newline);
        }

        private void WriteSections() {
            if (sections.Count == 0) {
                Section dummySection = new Section();
                dummySection.Level = 5;
                dummySection.ChartChangeOnField = this.ChartChangeOnField;
                dummySection.ChartLabelHeader = this.ChartLabelHeader;
                dummySection.ChartPercentageHeader = this.ChartPercentageHeader;
                dummySection.ChartShowAtBottom = this.ChartShowAtBottom;
                dummySection.ChartShowBorder = this.ChartShowAtBottom;
                dummySection.ChartTitle = this.ChartTitle;
                dummySection.ChartValueField = this.ChartValueField;
                dummySection.ChartValueHeader = this.ChartValueHeader;
                dummySection.IncludeChart = this.IncludeChart;
                htmlContent.Append("<TABLE Width='100%' class='TableStyle'  cellspacing=0 cellpadding=5 border=0>" + newline);
                if (this.IncludeChart && !this.ChartShowAtBottom)
                    GenerateBarChart("", dummySection);
                Hashtable total = WriteSectionDetail(null, "");
                if (this.IncludeTotal) {
                    dummySection.IncludeTotal = true;
                    WriteSectionFooter(dummySection, total);
                }
                if (this.IncludeChart && this.ChartShowAtBottom)
                    GenerateBarChart("", dummySection);
                htmlContent.Append("</TABLE></BODY></HTML>");
            }
            foreach (Section section in sections) {
                iLevel = 0;
                htmlContent.Append("<TABLE Width='100%' class='TableStyle'  cellspacing=0 cellpadding=5 border=0>" + newline);
                //RecurseSections(section, "");
                htmlContent.Append("</TABLE></BODY></HTML>");
            }
        }

        private void WriteSectionHeader(Section section, string sectionValue) {
            string bg = section.backColor;
            string style = " style=\"font-family: " + ReportFont + "; font-weight:bold; font-size:";
            style += getFontSize(section.Level);
            if (section.GradientBackground)
                style += "; " + gradientStyle.Replace("BackColor", bg) + "\"";
            else style += "\" bgcolor='" + bg + "' ";

            htmlContent.Append("<TR><TD colspan='" + this.ReportFields.Count + "' " + style + " >");
            htmlContent.Append(section.TitlePrefix + sectionValue);
            htmlContent.Append("</TD></TR>" + newline);
        }

        private Hashtable WriteSectionDetail(Section section, string criteria) {
            Hashtable totalArray = new Hashtable();
            totalArray = PrepareData(totalArray);
            if (section == null) {
                section = new Section();
            }
            try {
                //Draw DetailHeader
                htmlContent.Append("<TR>" + newline);
                string cellParams = "";
                foreach (Field field in this.reportFields) {
                    cellParams = " bgcolor='" + field.headerBackColor + "' ";
                    if (field.Width != 0)
                        cellParams += " WIDTH=" + field.Width + " ";
                    cellParams += " ALIGN='" + field.alignment + "' ";
                    htmlContent.Append("  <TD " + cellParams + " class='ColumnHeaderStyle'><b>" + field.HeaderName + "</b></TD>" + newline);
                }
                htmlContent.Append("</TR>" + newline);

                //Draw Data
                if (criteria == null || criteria.Trim() == "")
                    criteria = "";
                else
                    criteria = criteria.Substring(3);


                foreach (DataRow dr in reportSource.Tables[0].Select(criteria)) {
                    htmlContent.Append("<TR>" + newline);
                    foreach (Field field in this.reportFields) {
                        cellParams = " bgcolor='" + field.backColor + "' ";
                        if (field.Width != 0)
                            cellParams += " WIDTH=" + field.Width + " ";
                        //if total field, by default set to RIGHT align.
                        if (this.TotalFields.Contains(field.FieldName))
                            cellParams += " align='right' ";
                        cellParams += " ALIGN='" + field.alignment + "' ";
                        htmlContent.Append("  <TD " + cellParams + " VALIGN='top' class='DetailData'>" + dr[field.FieldName] + "</TD>" + newline);
                    }
                    htmlContent.Append("</TR>" + newline);
                    try {
                        foreach (object totalField in TotalFields) {
                            totalArray[totalField.ToString()] = float.Parse(totalArray[totalField.ToString()].ToString()) +
                                float.Parse(dr[totalField.ToString()].ToString());
                        }
                    }
                    catch (Exception exp) {
                        ;//to-do: show error message at total field
                    }
                }
            }
            catch (Exception err) {
                htmlContent.Append("<p align=CENTER><b>Unable to generate report.<br>" + err.Message + "</b></p>");
            }
            return totalArray;
        }

        private void WriteSectionFooter(Section section, Hashtable totalArray) {
            string cellParams = "";
            //Include Total row if specified.
            if (section.IncludeTotal) {
                htmlContent.Append("<TR>" + newline);
                foreach (Field field in this.reportFields) {
                    cellParams = "";
                    if (field.Width != 0)
                        cellParams += " WIDTH=" + field.Width + " ";
                    cellParams += " style=\"font-family: " + ReportFont + "; font-size:";
                    cellParams += getFontSize(section.Level + 1) + "; border-style:outline; border-width:1 \" ";
                    if (totalArray.Contains(field.FieldName)) {
                        htmlContent.Append("  <TD " + cellParams + " align='right'><u>Total: " + totalArray[field.FieldName].ToString() + "</u></TD> " + newline);
                    }
                    else {
                        htmlContent.Append("  <TD " + cellParams + ">&nbsp;</TD>" + newline);
                    }
                }
                htmlContent.Append("</TR>");
            }
        }

        private void WriteFooter() {
            htmlContent.Append("<BR>");
        }

        private Hashtable RecurseSections(Section section, string criteria) {
            iLevel++;
            section.Level = iLevel;
            ArrayList result = GetDistinctValues(this.reportSource, section.GroupBy, criteria);
            Hashtable ht = new Hashtable();
            PrepareData(ht);
            foreach (object obj in result) {
                Hashtable sectionTotal = new Hashtable();
                PrepareData(sectionTotal);
                //Construct critiera string to select data for the current section
                string tcriteria = criteria + "and " + section.GroupBy + "='" + obj.ToString().Replace("'", "''") + "' ";
                WriteSectionHeader(section, obj.ToString());
                //If user not specified to display chart at bottom of the section
                if (section.IncludeChart && !section.ChartShowAtBottom && !section.isChartCreated)
                    GenerateBarChart(tcriteria, section);
                if (section.SubSection != null) {
                    sectionTotal = RecurseSections(section.SubSection, tcriteria);
                    iLevel--;
                }
                else {
                    sectionTotal = WriteSectionDetail(section, tcriteria);
                    ht = AccumulateTotal(ht, sectionTotal);
                }
                //If user specified to display chart at bottom of the section
                WriteSectionFooter(section, sectionTotal);
                if (section.IncludeChart && section.ChartShowAtBottom && !section.isChartCreated)
                    GenerateBarChart(tcriteria, section);
                section.isChartCreated = false;
            }
            if (section.Level < 2)
                htmlContent.Append("<TR><TD colspan='" + this.ReportFields.Count + "'>&nbsp;</TD></TR>");

            return ht;
        }

        private void GenerateBarChart(string criteria, Section section) {
            string changeOnField = section.ChartChangeOnField;
            string valueField = section.ChartValueField;
            bool showBorder = section.ChartShowBorder;
            section.isChartCreated = true;
            string[] colors = { "#ff0000", "#ffff00", "#ff00ff", "#00ff00", "#00ffff", "#0000ff", "#ff0f0f", "#f0f000", "#ff00f0", "#0f00f0" };
            htmlContent.Append("<TR><TD colspan='" + this.ReportFields.Count + "' align=CENTER>" + newline);
            htmlContent.Append("<!--- Chart Table starts here -->" + newline);
            if (showBorder) {
                htmlContent.Append("<TABLE cellpadding=4 cellspacing=1 border=0 bgcolor='#f5f5f5' width=550> ");
            }
            else {
                htmlContent.Append("<TABLE border=0 cellspacing=5 width=550>");
            }
            if (criteria.ToUpper().StartsWith(" AND ")) {
                criteria = criteria.Substring(3);
            }
            try {
                ArrayList result = GetDistinctValuesForChart(this.reportSource, criteria, changeOnField, valueField);
                ArrayList labels = (ArrayList)result[0];
                ArrayList values = (ArrayList)result[1];
                float total = 0;
                foreach (Object obj in values) {
                    total += float.Parse(obj.ToString());
                }
                int ChartWidth = 300;

                string barTitle = "<TR><TD class='DetailHeader' colspan=3 align='CENTER' width=550><B>ChartTitle</B></TD></TR>";
                htmlContent.Append(barTitle.Replace("ChartTitle", section.ChartTitle) + newline);

                string barTemplate = "<TR><TD Width=150 class='DetailData' align='right'>Label</TD>" + newline;
                barTemplate += " <TD  class='DetailData' Width=" + (ChartWidth + 25) + ">" + newline;
                barTemplate += "    <TABLE cellpadding=0 cellspacing=0 HEIGHT='20' WIDTH=" + ChartWidth + " class='TableStyle'>" + newline;
                barTemplate += "       <TR>" + newline;
                barTemplate += "          <TD Width=ChartWidth>" + newline;
                barTemplate += "             <TABLE class='TableStyle' HEIGHT='20' Width=ChartTWidth border=NOTZERO>" + newline;
                barTemplate += "                <TR>" + newline;
                barTemplate += "                   <TD Width=ChartWidth bgcolor='BackColor' Width=ChartWidth style=\"FILTER: progid:DXImageTransform.Microsoft.Gradient(gradientType=0,startColorStr=BackColor,endColorStr=#ffffff); \"></TD>" + newline;
                barTemplate += "                </TR>" + newline;
                barTemplate += "             </TABLE>" + newline;
                barTemplate += "          </TD>" + newline;
                barTemplate += "          <TD class='DetailData'>&nbsp;Percentage</TD>" + newline;
                barTemplate += "       </TR>" + newline;
                barTemplate += "    </TABLE>";
                barTemplate += "</TD><TD Width=70 class='DetailData'>Value</TD></TR>";

                string barHTemplate = "<TR>" + newline;
                barHTemplate += " <TD Width=150  class='DetailData' align='right' bgColor='#e5e5e5'>Label</TD>" + newline;
                barHTemplate += " <TD  bgColor='#e5e5e5' class='DetailData' Width=" + (ChartWidth + 25) + ">";
                barHTemplate += "Percentage</TD>" + newline;
                barHTemplate += " <TD Width=25  class='DetailData' bgColor='#e5e5e5'>Value</TD></TR>";
                barHTemplate = barHTemplate.Replace("Label", section.ChartLabelHeader);
                barHTemplate = barHTemplate.Replace("Percentage", section.ChartPercentageHeader);
                barHTemplate = barHTemplate.Replace("Value", section.ChartValueHeader);
                htmlContent.Append(barHTemplate + newline);

                string temp;
                float width = 0;
                float val = 0;
                float percent = 0;
                int cntColor = 0;
                for (int i = 0; i < labels.Count; i++) {
                    temp = barTemplate;
                    val = float.Parse(values[i].ToString());
                    width = float.Parse(values[i].ToString()) * ChartWidth / total;
                    percent = float.Parse(values[i].ToString()) * 100 / total;

                    temp = temp.Replace("Label", labels[i].ToString());
                    if (percent == 0.0) {
                        temp = temp.Replace("BackColor", "#f5f5f5");
                        temp = temp.Replace("NOTZERO", "0");
                    }
                    else {
                        temp = temp.Replace("BackColor", colors[cntColor]);
                        temp = temp.Replace("NOTZERO", "1");
                    }
                    temp = temp.Replace("ChartTWidth", Decimal.Round((Decimal)(width + 2), 0).ToString());
                    temp = temp.Replace("ChartWidth", Decimal.Round((Decimal)width, 0).ToString());
                    temp = temp.Replace("Value", val.ToString());
                    temp = temp.Replace("Percentage", Decimal.Round((decimal)percent, 2).ToString() + "%");

                    htmlContent.Append(temp + newline);
                    cntColor++;
                    if (cntColor == 10)
                        cntColor = 0;
                }
            }
            catch (Exception err) {
                htmlContent.Append("<TR><TD valign=MIDDLE align=CENTER><b>Unable to Generate Chart.<br>" + err.Message + "</b></TD></TR>");
            }
            htmlContent.Append("</TABLE>" + newline);
            htmlContent.Append("<!--- Chart Table ends here -->");
            htmlContent.Append("</TD></TR>");
        }


        private ArrayList GetDistinctValuesForChart(DataSet dataSet, string criteria, string changeOnField, string valueField) {
            ArrayList result = new ArrayList();
            ArrayList distinctValues = new ArrayList();
            if (criteria == null || criteria.Trim() == "") {
                criteria = "";
            }
            else {
                criteria = criteria.Substring(3);
            }
            foreach (DataRow dr in dataSet.Tables[0].Select(criteria)) {
                if (!distinctValues.Contains(dr[changeOnField].ToString())) {
                    distinctValues.Add(dr[changeOnField].ToString());
                }
            }
            ArrayList totalValues = new ArrayList();
            if (criteria.Trim() != "")
                criteria += " and ";
            foreach (object obj in distinctValues) {
                DataRow[] rows = reportSource.Tables[0].Select(criteria + changeOnField + "='" + obj.ToString().Replace("'", "''") + "' ");
                if (valueField.Trim().ToUpper() == "COUNT") {
                    totalValues.Add(rows.Length.ToString());
                }
                else {
                    float sum = 0;
                    foreach (DataRow row in rows)
                        sum += float.Parse(row[valueField].ToString());
                    totalValues.Add(sum.ToString());
                }
            }
            result.Add(distinctValues);
            result.Add(totalValues);
            return result;
        }


        private ArrayList GetDistinctValues(DataSet dataSet, string columnName, string criteria) {
            ArrayList distinctValues = new ArrayList();
            if (criteria == null || criteria.Trim() == "") {
                criteria = "";
            }
            else {
                criteria = criteria.Substring(3);
            }
            foreach (DataRow dr in dataSet.Tables[0].Select(criteria)) {
                if (!distinctValues.Contains(dr[columnName].ToString())) {
                    distinctValues.Add(dr[columnName].ToString());
                }
            }
            return distinctValues;
        }


        private Hashtable PrepareData(Hashtable totalArray) {
            foreach (object obj in TotalFields) {
                if (!totalArray.Contains(obj.ToString())) {
                    totalArray.Add(obj.ToString(), 0.0F);
                }
            }
            return totalArray;
        }

        private Hashtable AccumulateTotal(Hashtable totalTable1, Hashtable totalTable2) {
            foreach (object totalField in TotalFields) {
                totalTable1[totalField.ToString()] = float.Parse(totalTable1[totalField.ToString()].ToString()) +
                    float.Parse(totalTable2[totalField.ToString()].ToString());
            }
            return totalTable1;
        }

        private string getFontSize(int level) {
            string fontSize = "";
            switch (level) {
                case 1:
                    fontSize = "14pt";
                    break;
                case 2:
                    fontSize = "12pt";
                    break;
                case 3:
                    fontSize = "10pt";
                    break;
                default:
                    fontSize = "9pt";
                    break;
            }
            return fontSize;
        }
        #endregion

    }

    #region Inner Classes

    internal class Section {
        internal string GroupBy;
        internal string TitlePrefix;
        internal bool IncludeFooter;
        internal bool GradientBackground;

        internal bool IncludeTotal;
        internal Section SubSection;
        /// <summary>
        /// HTML Color code as string
        /// </summary>
        internal string backColor;
        internal Color cBackColor;
        internal int Level;
        internal bool isChartCreated;

        internal bool IncludeChart;
        internal string ChartTitle;
        internal bool ChartShowAtBottom;
        internal string ChartChangeOnField;
        internal string ChartValueField = "Count";
        internal bool ChartShowBorder;
        internal string ChartLabelHeader = "Label";
        internal string ChartPercentageHeader = "Percentage";
        internal string ChartValueHeader = "Value";

        internal Color BackColor {
            set { backColor = Util.GetHTMLColorString(value); cBackColor = value; }
            get { return cBackColor; }
        }

        internal Section() {
            SubSection = null;
            BackColor = Color.FromArgb(240, 240, 240);
            ChartValueField = "Count";
            GradientBackground = false;
            ChartTitle = "&nbsp;";
        }

        internal Section(string groupBy, string titlePrefix) {
            GroupBy = groupBy;
            TitlePrefix = titlePrefix;
            SubSection = null;
            BackColor = Color.FromArgb(240, 240, 240);
            ChartValueField = "Count";
            GradientBackground = false;
            ChartTitle = "&nbsp;";
        }

        internal Section(string groupBy, string titlePrefix, Color bgcolor) {
            GroupBy = groupBy;
            TitlePrefix = titlePrefix;
            SubSection = null;
            BackColor = bgcolor;
            ChartValueField = "Count";
            GradientBackground = false;
            ChartTitle = "&nbsp;";
        }
    }

    internal class Field {
        internal string FieldName;
        internal string HeaderName;
        
        internal string backColor;
        internal Color cBackColor;
        
        internal string headerBackColor;
        internal Color cHeaderBackColor;
        internal int Width;
        internal bool isTotalField = false;
        internal string alignment = "LEFT";

        internal ALIGN Alignment {
            set {
                switch (value) {
                    case ALIGN.LEFT: alignment = "LEFT";  break;
                    case ALIGN.RIGHT: alignment = "RIGHT"; break;
                    case ALIGN.CENTER: alignment = "CENTER"; break;
                    default: alignment = "LEFT"; break;
                }
            }
            get {
                switch (alignment) {
                    case "LEFT": return ALIGN.LEFT;
                    case "RIGHT":  return ALIGN.RIGHT;
                    case "CENTER":  return ALIGN.CENTER;
                    default: return ALIGN.LEFT;
                }
            }
        }


        internal Color BackColor {
            set { backColor = Util.GetHTMLColorString(value); cBackColor = value; }
            get { return cBackColor; }
        }

        internal Color HeaderBackColor {
            set { headerBackColor = Util.GetHTMLColorString(value); cHeaderBackColor = value; }
            get { return cHeaderBackColor; }
        }

        internal Field() {
            FieldName = "";
            HeaderName = "Column Header";
            BackColor = Color.White;
            Width = 0;
            HeaderBackColor = Color.White;
        }

        internal Field(string fieldName, string headerName) {
            FieldName = fieldName;
            HeaderName = headerName;
            BackColor = Color.White;
            Width = 0;
            HeaderBackColor = Color.White;
        }

        internal Field(string fieldName, string headerName, int width) {
            FieldName = fieldName;
            HeaderName = headerName;
            BackColor = Color.White;
            Width = width;
            HeaderBackColor = Color.White;
        }

        internal Field(string fieldName, string headerName, int width, Color bgcolor) {
            FieldName = fieldName;
            HeaderName = headerName;
            BackColor = Color.White;
            Width = width;
            BackColor = bgcolor;
            HeaderBackColor = Color.White;
        }

        internal Field(string fieldName, string headerName, int width, ALIGN TextAlignment) {
            FieldName = fieldName;
            HeaderName = headerName;
            BackColor = Color.White;
            Width = width;
            BackColor = Color.White;
            HeaderBackColor = Color.White;
            Alignment = TextAlignment;
        }

        internal Field(string fieldName, string headerName, int width, Color bgcolor, Color headerBgColor) {
            FieldName = fieldName;
            HeaderName = headerName;
            Width = width;
            BackColor = bgcolor;
            HeaderBackColor = headerBgColor;
        }

        internal Field(string fieldName, string headerName, Color bgcolor, Color headerBgColor) {
            FieldName = fieldName;
            HeaderName = headerName;
            Width = 0;
            BackColor = bgcolor;
            HeaderBackColor = headerBgColor;
        }

        internal Field(string fieldName, string headerName, Color headerBgColor) {
            FieldName = fieldName;
            HeaderName = headerName;
            Width = 0;
            BackColor = Color.White;
            HeaderBackColor = headerBgColor;
        }
    }

    internal enum ALIGN {
        LEFT = 0,
        RIGHT,
        CENTER
    }

    internal class Util {
        internal static string GetHTMLColorString(Color color) {
            if (color.IsNamedColor)
                return color.Name;
            else
                return "#" + color.R.ToString("X2") + color.G.ToString("X2") + color.B.ToString("X2");
        }
    }

    #endregion
}
