using Accord.Math;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using System;
using System.Linq;
using System.Numerics;
using System.Windows;

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
            int plotIndex = 0;

            double[] squareWave = GenerateSquareWave(4);
            double[] carrierSignal = GenerateCarrierSignal(32);
            double[] amplitudeModulatedSignal = AmplitudeModulation(squareWave, carrierSignal);
            double[] frequencyModulatedSignal = FrequencyModulation(squareWave, 20);
            double[] phaseModulatedSignal = PhaseModulation(squareWave, 30);

            // Square Wave
            PlotGraph(squareWave, "Square Wave Modulation", plotIndex++);

            // Carrier Signal
            PlotGraph(carrierSignal, "Carrier Wave Modulation", plotIndex++);

            // Amplitude Modulation
            PlotGraph(amplitudeModulatedSignal, "Amplitude Modulation", plotIndex++);

            // Frequency Modulation
            PlotGraph(frequencyModulatedSignal, "Frequency Modulation", plotIndex++);

            // Phase Modulation
            PlotGraph(phaseModulatedSignal, "Phase Modulation", plotIndex++);

            // Spectrum of Amplitude Modulated Signal
            var amplitudeModulatedSpectrum = CalculateSpectrum(amplitudeModulatedSignal);
            PlotSpectrumGraph(amplitudeModulatedSpectrum.Item1, amplitudeModulatedSpectrum.Item2, "Amplitude Modulated Spectrum", plotIndex++);

            // Spectrum of Frequency Modulated Signal
            var frequencyModulatedSpectrum = CalculateSpectrum(frequencyModulatedSignal);
            PlotSpectrumGraph(frequencyModulatedSpectrum.Item1, frequencyModulatedSpectrum.Item2, "Frequency Modulated Spectrum", plotIndex++);

            // Spectrum of Phase Modulated Signal
            var phaseModulatedSpectrum = CalculateSpectrum(phaseModulatedSignal);
            PlotSpectrumGraph(phaseModulatedSpectrum.Item1, phaseModulatedSpectrum.Item2, "Phase Modulated Spectrum", plotIndex++);

            // Filtered Spectrum of Amplitude Modulated Signal
            var filteredAmplitudeModulatedSpectrum = FilterSpectrum(amplitudeModulatedSpectrum.Item1, amplitudeModulatedSpectrum.Item2, 25, 40);
            PlotFilteredSpectrumGraph(amplitudeModulatedSpectrum.Item1, filteredAmplitudeModulatedSpectrum, "Filtered Spectrum of Amplitude Modulated", plotIndex++);

            // Filtered Spectrum of Frequency Modulated Signal
            var filteredFrequencyModulatedSpectrum = FilterSpectrum(frequencyModulatedSpectrum.Item1, frequencyModulatedSpectrum.Item2, 20, 80);
            PlotFilteredSpectrumGraph(frequencyModulatedSpectrum.Item1, filteredFrequencyModulatedSpectrum, "Filtered Spectrum of Frequency Modulated", plotIndex++);

            // Filtered Spectrum of Phase Modulated Signal
            var filteredPhaseModulatedSpectrum = FilterSpectrum(phaseModulatedSpectrum.Item1, phaseModulatedSpectrum.Item2, 25, 35);
            PlotFilteredSpectrumGraph(phaseModulatedSpectrum.Item1, filteredPhaseModulatedSpectrum, "Filtered Spectrum of Phase Modulated", plotIndex++);


            // Synthesized and Envelope Signal
            var synthesizedSignal = SynthesizeSignal(filteredAmplitudeModulatedSpectrum);
            var envelope = CalculateEnvelope(synthesizedSignal);
            PlotSynthesizedSignalAndEnvelope(synthesizedSignal, envelope, "Synthesized and Envelope", plotIndex++);

            ShowPlots();
        }

        private double[] SynthesizeSignal(double[] filteredSpectrum)
        {
            Complex[] fftInput = filteredSpectrum.Select(x => new Complex(x, 0)).ToArray();
            FourierTransform.FFT(fftInput, FourierTransform.Direction.Backward);

            // Вещественная часть обратного преобразования Фурье
            return fftInput.Select(x => x.Real).ToArray();
        }

        private double[] CalculateEnvelope(double[] signal)
        {
            Complex[] complexSignal = signal.Select(x => new Complex(x, 0)).ToArray();

            // Метод FHT из HilbertTransform
            Accord.Math.HilbertTransform.FHT(complexSignal, FourierTransform.Direction.Forward);

            return complexSignal.Select(x => x.Magnitude).ToArray();
        }
        private (double[], double[]) CalculateSpectrum(double[] signal)
        {
            int n = signal.Length;
            int fftSize = 1;
            while (fftSize < n)
                fftSize *= 2;

            if (n < fftSize)
                Array.Resize(ref signal, fftSize);

            Complex[] fftInput = signal.Select(x => new Complex(x, 0)).ToArray();
            FourierTransform.FFT(fftInput, FourierTransform.Direction.Forward);

            double[] frequencies = new double[fftSize / 2];
            double[] magnitudes = new double[fftSize / 2];

            for (int i = 0; i < fftSize / 2; i++)
            {
                frequencies[i] = i * SamplingRate / fftSize;
                magnitudes[i] = 2 * fftInput[i].Magnitude / fftSize;
            }

            return (frequencies, magnitudes);
        }

        private void PlotSynthesizedSignalAndEnvelope(double[] synthesizedSignal, double[] envelope, string title, int plotIndex)
        {
            if (plotIndex >= 0 && plotIndex < plotModels.Length)
            {
                var signalSeries = new LineSeries();
                var envelopeSeries = new LineSeries();

                for (int i = 0; i < synthesizedSignal.Length; i++)
                {
                    signalSeries.Points.Add(new DataPoint(timeValues[i], synthesizedSignal[i]));
                    envelopeSeries.Points.Add(new DataPoint(timeValues[i], envelope[i]));
                }

                var signalModel = plotModels[plotIndex];
                signalModel.Series.Add(signalSeries);
                signalModel.Series.Add(envelopeSeries);
                signalModel.Title = title;
            }
            else
                Console.WriteLine($"Invalid plot index: {plotIndex}");
        }
        private void PlotGraph(double[] data, string title, int plotIndex)
        {
            if (plotIndex >= 0 && plotIndex < plotModels.Length)
            {
                var signalSeries = new LineSeries();
                for (int i = 0; i < data.Length; i++)
                    signalSeries.Points.Add(new DataPoint(timeValues[i], data[i]));

                var signalModel = plotModels[plotIndex];
                signalModel.Series.Add(signalSeries);
                signalModel.Title = title;
            }
            else
                Console.WriteLine($"Invalid plot index: {plotIndex}");
        }

        private void PlotSpectrumGraph(double[] frequencies, double[] magnitudes, string title, int plotIndex)
        {
            if (plotIndex >= 0 && plotIndex < plotModels.Length)
            {
                var spectrumSeries = new LineSeries();
                for (int i = 0; i < frequencies.Length; i++)
                    spectrumSeries.Points.Add(new DataPoint(frequencies[i], magnitudes[i]));

                var spectrumModel = plotModels[plotIndex];
                spectrumModel.Series.Add(spectrumSeries);
                spectrumModel.Title = title;

                // Настройка осей для центрирования
                spectrumModel.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Minimum = 0, Maximum = 100 });
                spectrumModel.Axes.Add(new LinearAxis { Position = AxisPosition.Left });
            }
            else
                Console.WriteLine($"Invalid spectrum plot index: {plotIndex}");
        }

        private void PlotFilteredSpectrumGraph(double[] frequencies, double[] filteredSpectrum, string title, int plotIndex)
        {
            if (plotIndex >= 0 && plotIndex < plotModels.Length)
            {
                var filteredSpectrumSeries = new LineSeries();
                for (int i = 0; i < frequencies.Length; i++)
                    filteredSpectrumSeries.Points.Add(new DataPoint(frequencies[i], filteredSpectrum[i]));

                var filteredSpectrumModel = plotModels[plotIndex];
                filteredSpectrumModel.Series.Add(filteredSpectrumSeries);
                filteredSpectrumModel.Title = title;
                filteredSpectrumModel.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "Frequency (Hz)", Minimum = 0, Maximum = 100 });
                filteredSpectrumModel.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Magnitude" });
            }
            else
                Console.WriteLine($"Invalid filtered spectrum plot index: {plotIndex}");
        }
        private void ShowPlots()
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
                result[i] = 0.5 * (Math.Sign(Math.Sin(2 * Math.PI * timeValues[i] * frequency)) + 1);

            return result;
        }

        private double[] GenerateCarrierSignal(double frequency)
        {
            double[] result = new double[timeValues.Length];
            for (int i = 0; i < timeValues.Length; i++)
                result[i] = Math.Sin(2 * Math.PI * frequency * timeValues[i]);

            return result;
        }

        private double[] AmplitudeModulation(double[] firstSignal, double[] secondSignal)
        {
            double[] result = new double[firstSignal.Length];
            for (int i = 0; i < firstSignal.Length; i++)
                result[i] = firstSignal[i] * secondSignal[i];
            return result;
        }

        private double[] FrequencyModulation(double[] signal, double modCoefficient)
        {
            double[] result = new double[timeValues.Length];
            for (int i = 0; i < timeValues.Length; i++)
                result[i] = Math.Sin(4 * Math.PI * (15 + modCoefficient * signal[i]) * timeValues[i]);
            return result;
        }

        private double[] PhaseModulation(double[] signal, double modCoefficient)
        {
            double[] result = new double[timeValues.Length];
            for (int i = 0; i < timeValues.Length; i++)
                result[i] = Math.Sin(2 * Math.PI * modCoefficient * timeValues[i] - 3 * modCoefficient * signal[i]);
            return result;
        }


        private double[] FilterSpectrum(double[] frequencies, double[] spectrum, double lowThreshold, double highThreshold)
        {
            // Пороговая фильтрация
            double[] filteredSpectrum = new double[spectrum.Length];

            for (int i = 0; i < spectrum.Length; i++)
            {
                if (Math.Abs(frequencies[i]) >= lowThreshold && Math.Abs(frequencies[i]) <= highThreshold)
                    filteredSpectrum[i] = spectrum[i];
                else
                    filteredSpectrum[i] = 0;
            }

            return filteredSpectrum;
        }
    }
}
