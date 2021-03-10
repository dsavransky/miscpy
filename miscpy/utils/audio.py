import numpy as np
import wavio
import warnings


class audio(object):
    """audio generation methods

    Args:
        samplingFrequency (int):
            Sampling Frequency in Hz (defaults to 44100).

    """

    def __init__(self, samplingFrequency=44100):

        assert isinstance(samplingFrequency, int)
        self.samplingFrequency = samplingFrequency

    def noteFreq(self, note, sharp=False, flat=False, octave=0):
        """
        Convert note name to frequency

        Args:
            note (char):
                Must be one of CDEFGAB
            sharp (bool):
                Set true for sharp
            flat (bool):
                Set true for flat
            ocatve (int):
                Octave specification. 0 represents octave starting with middle C.
                Negative values go lower, positive values go higher.  Octave -4
                represents the bottom of the standard 88 key piano such that
                ``noteFreq('A',octave=-4)`` is the lowest pitch (27.5 Hz).  Octave +3 is
                the top of a standard piano, such that ``noteFreq('C',octave=5)`` is the
                highest pitch (4186 Hz).
        Returns:
            float: Frequency in Hz

        """

        notenames = np.array(["C", "D", "E", "F", "G", "A", "B"])
        notenums = np.array([1, 3, 5, 6, 8, 10, 12])

        assert not (sharp and flat)
        assert note in notenames

        noteval = notenums[notenames == note][0] + 3 + (octave + 3) * 12
        if flat:
            noteval -= 1
        if sharp:
            noteval += 1

        return 2 ** ((noteval - 49) / 12) * 440

    def tone(self, freq, duration=0.5, e=0.1):
        """
        Generate sine wave corresponding to a frequency

        Args:
            freq (float):
                Frequency in Hz
            duration (float):
                Length of tone in seconds
            e (float):
                Window parameter.  If non-zero, a Planck-Tukey window will be applied
                to the tone.  Defults to 0.1


        Returns:
            numpy ndarray: The tone

        """

        t = np.linspace(
            0,
            duration,
            int(np.round(duration * self.samplingFrequency)),
            endpoint=False,
        )

        out = np.sin(2 * np.pi * freq * t)

        if e > 0:
            out *= self.planckTaperWindow(len(out), e=e)

        return out

    def writeWav(self, stream, filename):
        """
        Write 24 bit wav file

        Args:
            stream (numpy ndarray of floats):
                Content to write
            filename (str):
                Full path to write to

        Returns:
            None

        """

        wavio.write(filename, stream, self.samplingFrequency, sampwidth=3)

    def planckTaperWindow(self, N, e=0.1):
        """
        Calculate the Planck-Taper window for N samples

        Args:
            N (int):
                Length of window
            e (float):
                Shape paramter.  Must be between 0 and 1.  0 Returns tophat

        Returns:
            numpy ndarray of floats: Window function

        """
        out = np.zeros(N)
        ns = np.linspace(0, 1, N)
        inds = (ns < e) & (ns > 0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out[inds] = 1 / (1 + np.exp(e / ns[inds] - e / (e - ns[inds])))
        out[(ns >= e) & (ns < 0.5)] = 1
        tmp = out[ns < 0.5][::-1]
        if np.mod(N, 2) == 1:
            tmp = tmp[:-1]
        out[ns >= 0.5] = tmp

        return out
