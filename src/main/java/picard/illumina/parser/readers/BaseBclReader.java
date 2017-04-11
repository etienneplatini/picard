package picard.illumina.parser.readers;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.illumina.parser.BclData;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

class BaseBclReader {
    private static final byte BASE_MASK = 0x0003;
    private static final byte[] BASE_LOOKUP = new byte[]{'A', 'C', 'G', 'T'};
    final InputStream[] streams;
    final File[] streamFiles;
    final int[] outputLengths;
    private BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    final int[] numClustersPerCycle;
    final int cycles;

    BaseBclReader(int[] outputLengths, BclQualityEvaluationStrategy bclQualityEvaluationStrategy) {
        this.outputLengths = outputLengths;
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;

        int cycles = 0;
        for (final int outputLength : outputLengths) {
            cycles += outputLength;
        }

        this.cycles = cycles;
        this.streams = new InputStream[cycles];
        this.streamFiles = new File[cycles];
        this.numClustersPerCycle = new int[cycles];
    }

    BaseBclReader(int[] outputLengths) {
        this.outputLengths = outputLengths;

        int cycles = 0;
        for (final int outputLength : outputLengths) {
            cycles += outputLength;
        }

        this.cycles = cycles;
        this.streams = new InputStream[cycles];
        this.streamFiles = new File[cycles];
        this.numClustersPerCycle = new int[cycles];
    }

    int getNumCycles() {
        return cycles;
    }

    int getNumClusters() {
        return numClustersPerCycle[0];
    }

    InputStream open(final File file, final boolean seekable, final boolean isGzip, final boolean isBgzf) throws IOException {
        final String filePath = file.getAbsolutePath();

        try {
            // Open up a buffered stream to read from the file and optionally wrap it in a gzip stream
            // if necessary
            if (isBgzf) {
                // Only BlockCompressedInputStreams can seek, and only if they are fed a SeekableStream.
                return new BlockCompressedInputStream(IOUtil.maybeBufferedSeekableStream(file));
            } else if (isGzip) {
                if (seekable) {
                    throw new IllegalArgumentException(
                            String.format("Cannot create a seekable reader for gzip bcl: %s.", filePath)
                    );
                }
                return (IOUtil.maybeBufferInputStream(new GZIPInputStream(new FileInputStream(file), Defaults.BUFFER_SIZE / 2),
                        Defaults.BUFFER_SIZE / 2));
            } else {
                if (seekable) {
                    throw new IllegalArgumentException(
                            String.format("Cannot create a seekable reader for provided bcl: %s.", filePath)
                    );
                }
                return IOUtil.maybeBufferInputStream(new FileInputStream(file));
            }
        } catch (final FileNotFoundException fnfe) {
            throw new PicardException("File not found: (" + filePath + ")", fnfe);
        } catch (final IOException ioe) {
            throw new PicardException("Error reading file: (" + filePath + ")", ioe);
        }
    }

    void decodeBasecall(BclData bclData, int read, int cycle, int byteToDecode) {
        if (byteToDecode == 0) {
            bclData.bases[read][cycle] = (byte) '.';
            bclData.qualities[read][cycle] = (byte) 2;
        } else {
            bclData.bases[read][cycle] = BASE_LOOKUP[byteToDecode & BASE_MASK];
            bclData.qualities[read][cycle] = bclQualityEvaluationStrategy.reviseAndConditionallyLogQuality((byte) (byteToDecode >>> 2));
        }
    }

    void decodeQualityBinnedBasecall(BclData bclData, int read, int cycle, int byteToDecode, CycleData cycleData) {
        if (byteToDecode == 0) {
            bclData.bases[read][cycle] = (byte) '.';
            bclData.qualities[read][cycle] = 2;
        } else {
            bclData.bases[read][cycle] = BASE_LOOKUP[byteToDecode & BASE_MASK];
            bclData.qualities[read][cycle] = cycleData.qualityBins[byteToDecode >>> 2];
        }
    }

    class CycleData {
        final short version;
        final int headerSize;

        final byte bitsPerBasecall;
        final byte bitsPerQualityScore;
        final int numberOfBins;
        final byte[] qualityBins;

        final int numTiles;
        final TileData[] tileInfo;
        final boolean pfExcluded;

        CycleData(short version, int headerSize, byte bitsPerBasecall, byte bitsPerQualityScore, int numberOfBins,
                  byte[] qualityBins, int numTiles, TileData[] tileInfo, boolean pfExcluded) {
            this.version = version;
            this.headerSize = headerSize;
            this.bitsPerBasecall = bitsPerBasecall;
            this.bitsPerQualityScore = bitsPerQualityScore;
            this.numberOfBins = numberOfBins;
            this.qualityBins = qualityBins;
            this.numTiles = numTiles;
            this.tileInfo = tileInfo;
            this.pfExcluded = pfExcluded;
        }


    }


    class TileData {
        final int tileNum;
        final int numClustersInTile;
        final int uncompressedBlockSize;
        final int compressedBlockSize;

        TileData(int tileNum, int numClustersInTile, int uncompressedBlockSize, int compressedBlockSize) {
            this.tileNum = tileNum;
            this.numClustersInTile = numClustersInTile;
            this.uncompressedBlockSize = uncompressedBlockSize;
            this.compressedBlockSize = compressedBlockSize;
        }
    }
}
