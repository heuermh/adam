package org.bdgenomics.adam.serialization;

import com.esotericsoftware.kryo.Kryo;

class ADAMJavaKryoRegistrator extends org.disq_bio.disq.serializer.DisqKryoRegistrator {
    @Override
    public void registerClasses(final Kryo kryo) {
        super.registerClasses(kryo);
        kryo.register(org.apache.spark.sql.execution.datasources.InMemoryFileIndex.SerializableBlockLocation[].class);
    }
}
