using Accord.Math;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using System;
using System.Linq;
using System.Numerics;
using System.Windows;
using System.Windows.Controls.Primitives;

namespace SignalModeling
{
    public partial class MainWindow : Window
    {
        private const int SamplingRate = 1000;
        private const int TimePoints = 1000;

        private readonly double[] timeValues;
        private readonly PlotModel[] plotModels;

        public MainWindow()
        {
            InitializeComponent();
            timeValues = GenerateTimeArray();
            plotModels = new PlotModel[12];

            for (int i = 0; i < plotModels.Length; i++)
            {
                plotModels[i] = new PlotModel();
                plotModels[i].Title = $"Plot {i + 1}";
            }

            PlotSignals();
        }

        private void PlotSignals()
        {
         
            double[] modCoefficients = { 1, 20, 30 };

            int plotIndex = 0;

                double[] squareWave = GenerateSquareWave(4);
                double[] carrierSignal = GenerateCarrierSignal(32);

                // Square Wave Modeling
                PlotSignalAndSpectrum(squareWave, "Square Wave Modeling", plotIndex++);

                // Carrier Signal Modeling
                PlotSignalAndSpectrum(carrierSignal, "Carrier Signal Modeling", plotIndex++);

                foreach (var modCoefficient in modCoefficients)
                {
                    // Amplitude Modulation
                    double[] amplitudeModulatedSignal = AmplitudeModulation(squareWave, carrierSignal);
                    PlotSignalAndSpectrum(amplitudeModulatedSignal, "Amplitude Modulation", plotIndex++);

                    // Frequency Modulation
                    double[] frequencyModulatedSignal = FrequencyModulation(squareWave, modCoefficient);
                    PlotSignalAndSpectrum(frequencyModulatedSignal, "Frequency Modulation", plotIndex++);

                    // Phase Modulation
                    double[] phaseModulatedSignal = PhaseModulation(squareWave, modCoefficient);
                    PlotSignalAndSpectrum(phaseModulatedSignal, "Phase Modulation", plotIndex++);
                }
            
        }


        private void PlotSignalAndSpectrum(double[] signal, string modulationType, int plotIndex)
        {
            if (plotIndex < plotModels.Length && plotIndex + 6 < plotModels.Length)
            {
                var signalSeries = new LineSeries();
                for (int i = 0; i < timeValues.Length; i++)
                {
                    signalSeries.Points.Add(new DataPoint(timeValues[i], signal[i]));
                }

                var signalModel = plotModels[plotIndex];
                signalModel.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "Time (s)" });
                signalModel.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Amplitude" });
                signalModel.Series.Add(signalSeries);
                signalModel.Title = $"{modulationType} Signal";
                signalModel.Subtitle = $"Frequency: {modulationType.Contains("Frequency")}, Mod. Coeff.: {modulationType.Contains("Modulation")}";

                var spectrum = CalculateSpectrum(signal);
                var spectrumSeries = new LineSeries();
                for (int i = 0; i < spectrum.Item1.Length; i++)
                {
                    spectrumSeries.Points.Add(new DataPoint(spectrum.Item1[i], spectrum.Item2[i]));
                }

                var spectrumModel = plotModels[plotIndex + 6];
                spectrumModel.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "Frequency (Hz)", Minimum = 0, Maximum = 65 });
                spectrumModel.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Magnitude" });
                spectrumModel.Series.Add(spectrumSeries);
                spectrumModel.Title = $"{modulationType} Spectrum";
                spectrumModel.Subtitle = $"Frequency: {modulationType.Contains("Frequency")}, Mod. Coeff.: {modulationType.Contains("Modulation")}";
            }
            else
            {
                // Обработка сценария, когда индекс находится за пределами массива
                Console.WriteLine($"Invalid plot index: {plotIndex}");
            }
        }


        private double[] GenerateTimeArray()
        {
            double[] time = new double[TimePoints];
            for (int i = 0; i < TimePoints; i++)
                time[i] = (double)i / SamplingRate;

            return time;
        }

        private double[] GenerateSquareWave(double frequency)
        {
            double[] result = new double[timeValues.Length];
            for (int i = 0; i < timeValues.Length; i++)
            {
                result[i] = 0.5 * (Math.Sign(Math.Sin(2 * Math.PI * timeValues[i] * frequency)) + 1);
            }
            return result;
        }

        private double[] GenerateCarrierSignal(double frequency)
        {
            double[] result = new double[timeValues.Length];
            for (int i = 0; i < timeValues.Length; i++)
            {
                result[i] = Math.Sin(2 * Math.PI * frequency * timeValues[i]);
            }
            return result;
        }


        private double[] AmplitudeModulation(double[] sig1, double[] sig2)
        {
            // Поэлементное умножение массивов
            double[] result = new double[sig1.Length];
            for (int i = 0; i < sig1.Length; i++)
            {
                result[i] = sig1[i] * sig2[i];
            }
            return result;
        }


        private double[] FrequencyModulation(double[] signal, double modCoefficient)
        {
            double[] result = new double[timeValues.Length];
            for (int i = 0; i < timeValues.Length; i++)
            {
                result[i] = Math.Sin(4 * Math.PI * (15 + modCoefficient * signal[i]) * timeValues[i]);
            }
            return result;
        }



        private double[] PhaseModulation(double[] signal, double modCoefficient)
        {
            double[] result = new double[timeValues.Length];
            for (int i = 0; i < timeValues.Length; i++)
            {
                result[i] = Math.Sin(2 * Math.PI * modCoefficient * timeValues[i] - 3 * modCoefficient * signal[i]);
            }
            return result;
        }


        private (double[], double[]) CalculateSpectrum(double[] signal)
        {
            int n = signal.Length;
            int fftSize = 1;
            while (fftSize < n)
            {
                fftSize *= 2;
            }

            if (n < fftSize)
            {
                Array.Resize(ref signal, fftSize);
            }

            Complex[] fftInput = signal.Select(x => new Complex(x, 0)).ToArray();
            FourierTransform.FFT(fftInput, FourierTransform.Direction.Forward);

            double[] frequencies = new double[n / 2];
            for (int i = 0; i < n / 2; i++)
                frequencies[i] = i * SamplingRate / n;

            fftInput[0] = Complex.Zero;

            double[] magnitude = new double[fftSize / 2];
            for (int i = 0; i < fftSize / 2; i++)
                magnitude[i] = fftInput[i].Magnitude;

            return (frequencies, magnitude);
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            for (int i = 0; i < plotModels.Length; i++)
            {
                var plotView = new OxyPlot.Wpf.PlotView();
                plotView.Model = plotModels[i];
                plotView.Width = 430;
                plotView.Height = 253;

                uniformGrid.Children.Add(plotView);
            }
        }

        private void Window_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            MessageBox.Show("Closing event called");
            Closing += Window_Closing;

        }

    }
}