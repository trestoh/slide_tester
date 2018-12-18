﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace slide_tester
{
    class scan
    {
        public List<double> mz;
        public List<double> intensity;
    }

    public struct peak
    {
        public double mz;
        public double intensity;

        public peak(double m, double i)
        {
            mz = m;
            intensity = i;
        }

    }

    public struct subarrNode
    {
        public int front;
        public int back;
        public double score;

        public subarrNode(int f, int b, double s)
        {
            front = f;
            back = b;
            score = s;
        }

    }

    public struct roundedWindow
    {
        public double width;
        public double offset;
        public double upper;
        public double lower;

        public roundedWindow(double w, double o, double u, double l)
        {
            width = w;
            offset = o;
            upper = u;
            lower = l;
        }

    }

    class DynWinDriver
    {
        const double massDiff = 1.0033548378;

        public List<peak> get_target_range(List<peak> peaks, double mono_mz, int charge)
        {
            List<peak>.Enumerator copier = peaks.GetEnumerator();
            List<peak>.Enumerator look_ahead = peaks.GetEnumerator();


            List<peak> light_spec = new List<peak>();
            List<peak> iso_hits = new List<peak>();

            double[] target_isos = { mono_mz, (mono_mz + (massDiff / charge)), (mono_mz + 2 * (massDiff / charge)), (mono_mz + 3 * (massDiff / charge)), (mono_mz + 4 * (massDiff / charge)) };
            double tol = 20 * mono_mz * charge * (1 / 1000000.0);
            int curr_iso = 0;
            bool hit_iso = false;
            double iso_comp = target_isos[0];

            while (look_ahead.MoveNext() && look_ahead.Current.mz < (mono_mz - 0.1))
                copier.MoveNext();

            light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));

            while (copier.MoveNext() && copier.Current.mz < (mono_mz + 4 * (massDiff / charge) + 0.1))
            {
                if (curr_iso < 5 && Math.Abs(target_isos[curr_iso] - copier.Current.mz) < tol)
                {
                    iso_hits.Add(new peak(copier.Current.mz, copier.Current.intensity));
                    hit_iso = true;
                }

                else if (hit_iso)
                {
                    int best_match = 0;
                    double min_diff = Math.Abs(iso_hits[0].mz - target_isos[curr_iso]);
                    for (int i = 0; i < iso_hits.Count(); i++)
                    {
                        if (Math.Abs(iso_hits[i].mz - target_isos[curr_iso]) < min_diff)
                        {
                            best_match = i;
                            min_diff = Math.Abs(iso_hits[i].mz - target_isos[curr_iso]);
                        }
                    }

                    for (int i = 0; i < iso_hits.Count(); i++)
                    {
                        if (i == best_match)
                            light_spec.Add(new peak(iso_hits[i].mz, iso_hits[i].intensity));

                        else
                            light_spec.Add(new peak(iso_hits[i].mz, -1 * iso_hits[i].intensity));
                    }

                    iso_hits.Clear();
                    curr_iso++;

                    if (curr_iso == 5)
                    {
                        light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));
                        hit_iso = false;
                    }

                    else if (Math.Abs(target_isos[curr_iso] - copier.Current.mz) < tol)
                    {
                        iso_hits.Add(new peak(copier.Current.mz, copier.Current.intensity));
                        hit_iso = true;
                    }

                    else
                    {
                        light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));
                        hit_iso = false;
                    }

                }

                else
                {
                    light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));
                }

            }

            //new code to account for "full iso_hits" issue

            if (iso_hits.Count() > 0)
            {
                int best_match = 0;
                double min_diff = Math.Abs(iso_hits[0].mz - target_isos[curr_iso]);
                for (int i = 0; i < iso_hits.Count(); i++)
                {
                    if (Math.Abs(iso_hits[i].mz - target_isos[curr_iso]) < min_diff)
                    {
                        best_match = i;
                        min_diff = Math.Abs(iso_hits[i].mz - target_isos[curr_iso]);
                    }
                }

                for (int i = 0; i < iso_hits.Count(); i++)
                {
                    if (i == best_match)
                        light_spec.Add(new peak(iso_hits[i].mz, iso_hits[i].intensity));

                    else
                        light_spec.Add(new peak(iso_hits[i].mz, -1 * iso_hits[i].intensity));
                }
            }


            look_ahead = copier;
            if (look_ahead.MoveNext())
                light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));


            return light_spec;
        }

        public roundedWindow calculate_window(List<peak> light_spec, double mono_mz)
        {
            double temp_score;
            double min_window = 0.4;
            double default_width = 4.0;
            double margin = 0.1;
            int precision = 1;

            List<subarrNode> first_scores = new List<subarrNode>();
            first_scores.Add(new subarrNode(0, 0, light_spec[0].intensity));
            for (int i = 1; i < light_spec.Count; i++)
            {
                double new_score = light_spec[i].intensity + first_scores[i - 1].score;
                if (new_score > light_spec[i].intensity)
                    first_scores.Add(new subarrNode(first_scores[i - 1].front, i, new_score));

                else
                    first_scores.Add(new subarrNode(i, i, light_spec[i].intensity));

            }


            double local_sum = light_spec[0].intensity;
            int j = 0;

            subarrNode best = new subarrNode(0, 0, 0);

            while (j < light_spec.Count)
            {
                j++;
                if (light_spec[j].mz - light_spec[0].mz >= min_window)
                {
                    best.back = j;
                    break;
                }
                local_sum += light_spec[j].intensity;
            }

            best.score = local_sum;

            int front = best.front;
            int back = best.back;

            while (back < (light_spec.Count - 1))
            {
                double next_mz = light_spec[back + 1].mz;
                int num_peaks = (back - front) + 1;
                if (num_peaks > 1 && (next_mz - light_spec[front].mz) > min_window)
                {
                    local_sum -= light_spec[front].intensity;
                    front++;
                }

                else if ((next_mz - light_spec[front].mz) < min_window)
                {
                    local_sum += light_spec[back + 1].intensity;
                    back++;
                }

                else if ((next_mz - light_spec[front].mz) > min_window || (next_mz - light_spec[front].mz) == min_window)
                {
                    local_sum = local_sum + light_spec[back + 1].intensity - light_spec[front].intensity;
                    front++;
                    back++;
                }

                if (local_sum > best.score)
                    best = new subarrNode(front, back, local_sum);
                if (local_sum + first_scores[front - 1].score > best.score)
                    best = new subarrNode(first_scores[front - 1].front, back, local_sum + first_scores[front - 1].score);


            }


            double upper, lower;

            lower = Math.Min(light_spec[best.front - 1].mz + margin, (light_spec[best.front].mz + light_spec[best.front - 1].mz) / 2.0);
            upper = Math.Max(light_spec[best.back + 1].mz - margin, (light_spec[best.back].mz + light_spec[best.back + 1].mz) / 2.0);


            double offset = ((upper + lower) / 2.0) - mono_mz;
            double width = upper - lower;

            roundedWindow final_window = new roundedWindow(width, offset, upper, lower);
            return final_window;
        }


        public roundedWindow fix_window(roundedWindow orig_window, double mono_mz, int charge)
        {

            double offset = orig_window.offset;
            double width = orig_window.width;
            double upper = orig_window.upper;
            double lower = orig_window.lower;

            if (width > 4.0)
            {
                if (offset < 0.0)
                {
                    lower = upper - 4.0;
                }
                else if (offset > (4 * (massDiff / charge)))
                {
                    upper = lower + 4.0;
                }
                offset = ((upper + lower) / 2.0) - mono_mz;
                width = upper - lower;
            }

            roundedWindow final_window = new roundedWindow(width, offset, upper, lower);
            return final_window;
        }


    }


    class Program
    {
        const double massDiff = 1.0033548378;

        

        static void Main(string[] args)
        {
            string line;

            var driver = new DynWinDriver();

            System.IO.StreamReader file = new System.IO.StreamReader(@"C:\dev\real_time_seq\diw_testing\diw_testing\tests\standard_noisy.txt");
            List<peak> temp_scan = new List<peak>();

            while ((line = file.ReadLine()) != null)
            {
                if (line == "BREAK")
                    break;
                string[] peak = line.Split('\t');
                temp_scan.Add(new peak(System.Convert.ToDouble(peak[0]), System.Convert.ToDouble(peak[1])));

            }

            List<peak>.Enumerator copier = temp_scan.GetEnumerator();
            List<peak>.Enumerator look_ahead = temp_scan.GetEnumerator();

            List<peak> light_spec = new List<peak>();
            List<peak> iso_hits = new List<peak>();

            int charge = 3;
            string temp;
            double mono_mz = 528.930236816406;

            double[] target_isos = { mono_mz, (mono_mz + (massDiff / charge)), (mono_mz + 2 * (massDiff / charge)), (mono_mz + 3 * (massDiff / charge)), (mono_mz + 4 * (massDiff / charge)) };
            double tol = 20 * mono_mz * charge * (1 / 1000000.0);
            int curr_iso = 0;
            bool hit_iso = false;
            double iso_comp = target_isos[0];

            while (look_ahead.MoveNext() && look_ahead.Current.mz < (mono_mz - 0.1))
                copier.MoveNext();

            light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));

            while (copier.MoveNext() && copier.Current.mz < (mono_mz + 4 * (massDiff / charge) + 0.1))
            {
                if (curr_iso < 5 && Math.Abs(target_isos[curr_iso] - copier.Current.mz) < tol)
                {
                    iso_hits.Add(new peak(copier.Current.mz, copier.Current.intensity));
                    hit_iso = true;
                }

                else if (hit_iso)
                {
                    int best_match = 0;
                    double min_diff = Math.Abs(iso_hits[0].mz - target_isos[curr_iso]);
                    for (int i = 0; i < iso_hits.Count(); i++)
                    {
                        if (Math.Abs(iso_hits[i].mz - target_isos[curr_iso]) < min_diff)
                        {
                            best_match = i;
                            min_diff = Math.Abs(iso_hits[i].mz - target_isos[curr_iso]);
                        }
                    }

                    for (int i = 0; i < iso_hits.Count(); i++)
                    {
                        if (i == best_match)
                            light_spec.Add(new peak(iso_hits[i].mz, iso_hits[i].intensity));

                        else
                            light_spec.Add(new peak(iso_hits[i].mz, -1 * iso_hits[i].intensity));
                    }

                    iso_hits.Clear();
                    curr_iso++;

                    if (curr_iso == 5)
                    {
                        light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));
                        hit_iso = false;
                    }

                    else if (Math.Abs(target_isos[curr_iso] - copier.Current.mz) < tol)
                    {
                        iso_hits.Add(new peak(copier.Current.mz, copier.Current.intensity));
                        hit_iso = true;
                    }

                    else
                    {
                        light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));
                        hit_iso = false;
                    }

                }

                else
                {
                    light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));
                }

            }

            //new code to account for "full iso_hits" issue
            
            if(iso_hits.Count() > 0)
            {
                int best_match = 0;
                double min_diff = Math.Abs(iso_hits[0].mz - target_isos[curr_iso]);
                for (int i = 0; i < iso_hits.Count(); i++)
                {
                    if (Math.Abs(iso_hits[i].mz - target_isos[curr_iso]) < min_diff)
                    {
                        best_match = i;
                        min_diff = Math.Abs(iso_hits[i].mz - target_isos[curr_iso]);
                    }
                }

                for (int i = 0; i < iso_hits.Count(); i++)
                {
                    if (i == best_match)
                        light_spec.Add(new peak(iso_hits[i].mz, iso_hits[i].intensity));

                    else
                        light_spec.Add(new peak(iso_hits[i].mz, -1 * iso_hits[i].intensity));
                }
            }
            

            look_ahead = copier;
            if (look_ahead.MoveNext())
                light_spec.Add(new peak(copier.Current.mz, -1 * copier.Current.intensity));

            int front_index = 0;
            //IEnumerator<ICentroid> walker = temp_ms1.Centroids.GetEnumerator();
            //walker.MoveNext();

            //this is likely unnecessary now
            
            while (light_spec[front_index].mz < (mono_mz - 0.1))
            {
                front_index++;
            }

            List<peak> light_spec_2 = driver.get_target_range(temp_scan, mono_mz, charge);


            for (int i = 0; i < light_spec_2.Count; i++)
            {
                if ( Math.Abs(light_spec[i].mz - light_spec_2[i].mz) >= 0.0001 || Math.Abs(light_spec[i].intensity - light_spec_2[i].intensity) >= 0.0001)
                {
                    Console.WriteLine("Something is wrong with modular function");
                    System.Environment.Exit(1);
                }
            }


            int back_index = front_index;
            double temp_score;
            double min_window = 0.4;
            double default_width = 4.0;
            double margin = 0.1;
            int precision = 1;

            //old implementation
            /*
            double window_limit;
            
            if ((light_spec[front_index].mz - light_spec[front_index - 1].mz) > 0.2)
                window_limit = (light_spec[front_index].mz + light_spec[front_index - 1].mz) / 2.0;

            else
                window_limit = light_spec[front_index].mz - 0.1;

            while (back_index < light_spec.Count() && (light_spec[back_index].mz - window_limit) < 0.4)
                back_index++;

            int min_index = back_index;
            List<subarrNode> first_scores = new List<subarrNode>();
            while (back_index < (light_spec.Count()))
            {
                int temp_index = front_index;
                temp_score = 0.0;
                while (temp_index <= back_index)
                {
                    temp_score += light_spec[temp_index].intensity;
                    temp_index++;
                }

                first_scores.Add(new subarrNode(front_index, back_index, temp_score));
                back_index++;
                if (back_index >= light_spec.Count())
                    break;
                while ((light_spec[back_index].mz - light_spec[front_index].mz) >= 0.4)
                    front_index++;

                front_index--;
            }

            List<subarrNode> final_scores = new List<subarrNode>();
            back_index = min_index + 1;
            final_scores.Add(first_scores[0]);

            while (back_index < (light_spec.Count() - 1))
            {
                subarrNode previous = final_scores[back_index - min_index - 1];
                double extended = previous.score + light_spec[back_index].intensity;
                if (extended > first_scores[back_index - min_index].score)
                {
                    final_scores.Add(new subarrNode(previous.front, previous.back + 1, extended));
                }

                else
                    final_scores.Add(first_scores[back_index - min_index]);

                back_index++;

            }

            subarrNode best = final_scores[0];
            for (int i = 1; i < final_scores.Count(); i++)
            {
                if (final_scores[i].score > best.score)
                    best = final_scores[i];
            }
            */

            //NEW IMPLEMENTATION
            
            List<subarrNode> first_scores = new List<subarrNode>();
            first_scores.Add(new subarrNode(0, 0, light_spec[0].intensity));
            for (int i = 1; i < light_spec.Count; i++)
            {
                double new_score = light_spec[i].intensity + first_scores[i - 1].score;
                if (new_score > light_spec[i].intensity)
                    first_scores.Add(new subarrNode(first_scores[i - 1].front, i, new_score));

                else
                    first_scores.Add(new subarrNode(i, i, light_spec[i].intensity));

            }


            double local_sum = light_spec[0].intensity;
            int j = 0;

            subarrNode best = new subarrNode(0, 0, 0);

            while (j < light_spec.Count)
            {
                j++;
                if (light_spec[j].mz - light_spec[0].mz >= min_window)
                {
                    best.back = j;
                    break;
                }
                local_sum += light_spec[j].intensity;
            }

            best.score = local_sum; 

            int front = best.front;
            int back = best.back;



            while (back < (light_spec.Count - 1))
            {
                double next_mz = light_spec[back + 1].mz;
                int num_peaks = (back - front) + 1;
                if (num_peaks > 1 && (next_mz - light_spec[front].mz) > min_window)
                {
                    local_sum -= light_spec[front].intensity;
                    front++;
                }

                else if ((next_mz - light_spec[front].mz) < min_window)
                {
                    local_sum += light_spec[back + 1].intensity;
                    back++;
                }

                else if ((next_mz - light_spec[front].mz) > min_window || (next_mz - light_spec[front].mz) == min_window)
                {
                    local_sum = local_sum + light_spec[back + 1].intensity - light_spec[front].intensity;
                    front++;
                    back++;
                }

                if (local_sum > best.score)
                    best = new subarrNode(front, back, local_sum);
                if (local_sum + first_scores[front - 1].score > best.score)
                    best = new subarrNode(first_scores[front - 1].front, back, local_sum + first_scores[front - 1].score);


            }



            using (StreamWriter writetext = new StreamWriter("C:\\dev\\real_time_seq\\diw_testing\\diw_testing\\tests\\standard_noisy_answers.txt"))
            {
                
                for (int i = 0; i < light_spec.Count; i++)
                {
                    string spec = (light_spec[i].mz).ToString() + '\t' + (light_spec[i].intensity).ToString();
                    writetext.WriteLine(spec);
                }

                string break_spec = "BREAK";
                writetext.WriteLine(break_spec);
            }

            double upper, lower;

            /*
            if ((light_spec[best.front].mz - light_spec[best.front - 1].mz) > 0.2)
                lower = (light_spec[best.front].mz + light_spec[best.front - 1].mz) / 2.0;

            else
                lower = light_spec[best.front].mz - 0.1;

            if ((light_spec[best.back + 1].mz - light_spec[best.back].mz) > 0.2)
                upper = (light_spec[best.back].mz + light_spec[best.back + 1].mz) / 2.0;

            else
                upper = light_spec[best.back].mz + 0.1;

            if (upper - lower > 4.0)
                upper = lower + 4.0;

            double offset = ((upper + lower) / 2.0) - mono_mz;
            double width = upper - lower;
            */

            
            lower = Math.Min(light_spec[best.front - 1].mz + margin, (light_spec[best.front].mz + light_spec[best.front - 1].mz) / 2.0);
            
            upper = Math.Max(light_spec[best.back + 1].mz - margin, (light_spec[best.back].mz + light_spec[best.back + 1].mz) / 2.0);


            //this line should ensure windows no larger than 4.0 MZ
            //if (upper - lower > 4.0)
            //    upper = lower + 4.0;

            double offset = ((upper + lower) / 2.0) - mono_mz;
            double width = upper - lower;

            roundedWindow dyn_win = driver.calculate_window(light_spec_2, mono_mz);

            using (StreamWriter writetext = new StreamWriter("C:\\dev\\real_time_seq\\diw_testing\\diw_testing\\tests\\standard_noisy_answers.txt", true))
            {
                string window_rep = dyn_win.offset.ToString() + '\t' + dyn_win.width.ToString() + '\t' + dyn_win.upper.ToString() + '\t' + dyn_win.lower.ToString();
                writetext.WriteLine(window_rep);
            }

            if (Math.Abs(dyn_win.width - width) >= 0.0001 || Math.Abs(dyn_win.offset - offset) >= 0.0001)
            {
                Console.WriteLine("Something is wrong with modular function");
                System.Environment.Exit(1);
            }

            if (width > 4.0)
            {
                /*
                if (lower < target_isos[0] && upper > target_isos[4])
                {
                    lower = target_isos[2] - 2.0;
                    upper = target_isos[2] + 2.0;
                }
                    */
                /*else*/if (offset < 0.0)
                {
                    lower = upper - 4.0;
                }
                else if (offset > (4 * (massDiff / charge)))
                {
                    upper = lower + 4.0;
                }
                offset = ((upper + lower) / 2.0) - mono_mz;
                width = upper - lower;
            }


            dyn_win = driver.fix_window(dyn_win, mono_mz, charge);


            using (StreamWriter writetext = new StreamWriter("C:\\dev\\real_time_seq\\diw_testing\\diw_testing\\tests\\standard_noisy_answers.txt", true))
            {
                string window_rep = dyn_win.offset.ToString() + '\t' + dyn_win.width.ToString() + '\t' + dyn_win.upper.ToString() + '\t' + dyn_win.lower.ToString();
                writetext.WriteLine(window_rep);
            }


            //offset = 10.12356456;
            offset = -1.23856;

            Console.WriteLine("Original Offset: {0}", offset);

            double off_rounded = (double)(Math.Round((double)offset, 1));
            double width_rounded = (double)(Math.Round((double)width, 1));

            double less_rounded = (double)(Math.Round((double)offset, 5));
            double first_few = (double)(Math.Round((double)offset, 2));

            double low_mass = 200.0 + first_few * 10;
            double remain_offset = less_rounded - first_few;
            double high_mass = 2000.0 + remain_offset * 10000;

            double recon_offset = (high_mass - 2000.0) / 10000 + (low_mass - 200.0) / 10;

            Console.WriteLine("Low: {0}, High: {1}, Reconstructed Offset: {2}", low_mass, high_mass, recon_offset);
            Console.WriteLine("Offset: {0}, Width {1}", off_rounded, width_rounded);

            temp_score = 0.0;

        }
    }
}
