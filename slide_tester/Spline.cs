using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using System.Xml;
using System.Xml.Schema;

namespace Spline
{
    class IsotopeSplineXMLParser
    {

        //Base64.Decoder decoder;
        static XmlReader reader;

        public IsotopeSplineXMLParser()
        {
            //decoder = Base64.getDecoder();
        }

        public Dictionary<int, List<CubicSpline>> parse(String path)
        {

            try
            {
                reader = XmlReader.Create(path);

                XmlDocument doc = new XmlDocument();
                doc.Load(reader);


                Dictionary<int, List<CubicSpline>> models = parseDocument(doc);
                return models;
            }

            catch (Exception ex)
            {

            }

            return null;
        }

        private Dictionary<int, List<CubicSpline>> parseDocument(XmlDocument doc)
        {
            // get root element
            XmlElement root = doc.DocumentElement;

            //Element rootEle = dom.getDocumentElement();
            int numModels = System.Convert.ToInt32(root.GetAttribute("maxIsotopeDepth"));

            Dictionary<int, List<CubicSpline>> numSulfur2models = new Dictionary<int, List<CubicSpline>>();

            foreach (XmlElement node in root.GetElementsByTagName("model"))
            {


                int numSulfur = -1, isotope = 0;

                isotope = System.Convert.ToInt32(node.GetAttribute("isotope"));

                if (node.HasAttribute("S"))
                    numSulfur = System.Convert.ToInt32(node.GetAttribute("S"));


                if (!numSulfur2models.ContainsKey(numSulfur))
                {
                    numSulfur2models[numSulfur] = new List<CubicSpline>(numModels);
                    for (int i = 0; i < numModels; i++)
                        numSulfur2models[numSulfur].Add(null);
                }
                //Console.Write(numSulfur2models[numSulfur].Count);

                numSulfur2models[numSulfur][isotope] = parseModel(node);

            }

            return numSulfur2models;
        }

        private CubicSpline parseModel(XmlElement modelEle)
        {

            double[] knots = decodeDoubleList(getTextValue(modelEle, "knots"));
            double[] coefficients = decodeDoubleList(getTextValue(modelEle, "coefficients"));

            List<Double> a = new List<Double>();
            List<Double> b = new List<Double>();
            List<Double> c = new List<Double>();
            List<Double> d = new List<Double>();
            List<Double> x = new List<Double>();

            for (int i = 0; i < coefficients.GetLength(0); i += 4)
            {
                a.Add(coefficients[i]);
                b.Add(coefficients[i + 1]);
                c.Add(coefficients[i + 2]);
                d.Add(coefficients[i + 3]);
            }

            for (int i = 0; i < knots.GetLength(0); ++i)
            {
                x.Add(knots[i]);
            }

            return new CubicSpline(a, b, c, d, x);
        }

        private double[] decodeDoubleList(String encoded)
        {

            //byte[] AsBytes = new byte[encoded.Length];
            //Buffer.BlockCopy(encoded.ToCharArray(), 0, AsBytes, 0, AsBytes.Length);

            Byte[] bytes = Convert.FromBase64String(encoded);

            //ByteBuffer buff;
            //byte[] decoded = decoder.decode(encoded.getBytes());

            double[] values = new double[bytes.Length / 8];

            for (int i = 0; i < bytes.Length; i += 8)
            {
                byte[] data = new byte[8];
                for (int j = 0; j < 8; ++j)
                {
                    data[j] = bytes[i + j];
                }

                //buff = ByteBuffer.wrap(data);
                //buff = buff.order(ByteOrder.LITTLE_ENDIAN);

                //values[i / 8] = buff.getDouble();
                values[i / 8] = BitConverter.ToDouble(data, 0);

            }

            return values;
        }

        private String getTextValue(XmlElement ele, String tagName)
        {
            String textVal = null;
            XmlNodeList nl = ele.GetElementsByTagName(tagName);
            if (nl != null && nl.Count > 0)
            {
                XmlElement el = (XmlElement)nl[0];
                textVal = el.ChildNodes[0].Value;
            }

            return textVal;
        }


    }

    public class IsotopeDistribution
    {

        private List<Double> massData;
        private List<Double> intensityData;

        public IsotopeDistribution()
        {
            massData = new List<Double>();
            intensityData = new List<Double>();
        }

        public IsotopeDistribution(int maxIsotope)
        {
            massData = new List<Double>(maxIsotope + 1);
            intensityData = new List<Double>(maxIsotope + 1);
            for (int i = 0; i <= maxIsotope; ++i)
            {
                massData.Add(0.0);
                intensityData.Add(0.0);
            }
        }

        public int size()
        {
            return massData.Count;
        }

        public double getIntensity(int i)
        {
            return intensityData.ElementAt<Double>(i);
        }

        public void setIntensity(int i, double intensity)
        {
            intensityData[i] = intensity;
        }

        public void setMass(int i, double mz)
        {
            massData[i] = mz;
        }

        public double getMass(int i)
        {
            return massData.ElementAt<Double>(i);
        }

        public List<Double> getMassData()
        {
            return massData;
        }

        public List<Double> getIntensityData()
        {
            return intensityData;
        }

        public void normalizeToBasePeak()
        {
            normalizeToValue(intensityData.Max());
        }

        public void normalizeToValue(double value)
        {
            for (int i = 0; i < intensityData.Count; ++i)
            {
                intensityData[i] = (intensityData.ElementAt<Double>(i) / value);
            }
        }
    }

    class CubicSpline
    {

        private List<Double> a_; // constant spline coefficients
        private List<Double> b_; // linear spline coefficients
        private List<Double> c_; // quadratic spline coefficients
        private List<Double> d_; // cubic spline coefficients
        private List<Double> x_; // knots

        public CubicSpline()
        {
            a_ = new List<Double>();
            b_ = new List<Double>();
            c_ = new List<Double>();
            d_ = new List<Double>();
            x_ = new List<Double>();
        }

        public CubicSpline(List<Double> a, List<Double> b, List<Double> c, List<Double> d, List<Double> x)
        {
            a_ = a;
            b_ = b;
            c_ = c;
            d_ = d;
            x_ = x;
        }

        public double eval(double x)
        {

            //int index = Collections.binarySearch(x_, x);
            int index = x_.BinarySearch(x);

            if (index < 0)
            {
                index = -(index + 2);
            }

            double xx = x - x_.ElementAt<Double>(index);

            return ((d_.ElementAt<Double>(index) * xx + c_.ElementAt<Double>(index)) * xx + b_.ElementAt<Double>(index)) * xx + a_.ElementAt<Double>(index);
        }

        public bool inBounds(double x)
        {
            return x >= x_.ElementAt<Double>(0) && x <= x_.ElementAt<Double>(x_.Count - 1);
        }

    }

    class IsotopeSplineDB
    {

        private static double C13C12_MASSDIFF_U = 1.0033548378;

        // a key of -1 is for the average splines
        // the index in List<CublicSpline> is the isotope
        public Dictionary<int, List<CubicSpline>> numSulfur2models_;

        public IsotopeSplineDB(String splinePath)
        {
            readSplinesFromFile_(splinePath);
        }

        private void readSplinesFromFile_(String splinePath)
        {
            IsotopeSplineXMLParser parser = new IsotopeSplineXMLParser();

            numSulfur2models_ = parser.parse(splinePath);
        }

        public IsotopeDistribution estimateFromPeptideWeight(double monoMass, int maxIsotope)
        {
            return estimateFromPeptideWeightAndSulfur(monoMass, maxIsotope, -1);
        }

        IsotopeDistribution estimateFromPeptideWeightAndSulfur(double monoMass, int maxIsotope, int numSulfur)
        {
            IsotopeDistribution result;

            if (inModelBounds(monoMass, maxIsotope, numSulfur))
            {
                result = new IsotopeDistribution(maxIsotope);
                for (int isotope = 0; isotope <= maxIsotope; ++isotope)
                {
                    //double probability = numSulfur2models_.get(numSulfur).get(isotope).eval(monoMass);
                    double probability = numSulfur2models_[numSulfur][isotope].eval(monoMass);
                    result.setMass(isotope, monoMass + (isotope * C13C12_MASSDIFF_U));
                    result.setIntensity(isotope, (float)probability);
                }
            }
            else
            {
                throw new Exception("Request out of range.");
            }

            return result;
        }

        public IsotopeDistribution estimateForFragmentFromWeights(double monoPeptideMass, double monoFragmentMass,
                                                                  int minIsotope, int maxIsotope)
        {

            return estimateForFragmentFromWeightsAndSulfur(monoPeptideMass, monoFragmentMass, minIsotope, maxIsotope, -1, -1);
        }

        public IsotopeDistribution estimateForFragmentFromWeightsAndSulfur(double monoPeptideMass, double monoFragmentMass,
                                                                  int minIsotope, int maxIsotope, int precursorSulfur, int fragmentSulfur)
        {

            IsotopeDistribution fragment = estimateFromPeptideWeightAndSulfur(monoFragmentMass, maxIsotope, precursorSulfur);
            IsotopeDistribution compFragment = estimateFromPeptideWeightAndSulfur(monoPeptideMass - monoFragmentMass, maxIsotope, fragmentSulfur);
            return calcFragmentIsotopeDistribution(fragment, compFragment, minIsotope, maxIsotope);

        }

        public bool inModelBounds(double monoMass, int maxIsotope, int numSulfur)
        {
            // Check if we have a sulfur-specific model for this
            if (!numSulfur2models_.ContainsKey(numSulfur))
            {
                return false;
            }


            // Check if max isotope is in bounds
            if (maxIsotope > numSulfur2models_[numSulfur].Count - 1)
            {
                return false;
            }

            // Check if masses are in bounds
            for (int isotope = 0; isotope <= maxIsotope; ++isotope)
            {
                if (!numSulfur2models_[numSulfur][isotope].inBounds(monoMass))
                {
                    return false;
                }
            }

            // All checks passed
            return true;
        }

        private IsotopeDistribution calcFragmentIsotopeDistribution(IsotopeDistribution fragment, IsotopeDistribution compFragment,
                                                                    int minIsotope, int maxIsotope)
        {

            IsotopeDistribution result = new IsotopeDistribution(maxIsotope);

            for (int i = 0; i < fragment.size(); ++i)
            {
                for (int isotope = minIsotope; isotope <= maxIsotope; ++isotope)
                {
                    if (isotope >= i && (isotope - i) < compFragment.size())
                    {
                        result.setIntensity(i, (float)(result.getIntensity(i) + compFragment.getIntensity(isotope - i)));
                    }
                }
                result.setIntensity(i, (float)(result.getIntensity(i) * fragment.getIntensity(i)));
                result.setMass(i, fragment.getMass(0) + (i * C13C12_MASSDIFF_U));
            }

            return result;
        }

    }
}
