import test from 'ava'
import path from 'path'

import { readImageNode, writeImageNode } from '@itk-wasm/image-io'

import { morphologicalContourInterpolationNode } from '../../dist/index-node.js'
import { testInputPath, testOutputPath } from './common.js'

test('Test morphologicalContourInterpolationNode', async t => {
  const testInputFilePath = path.join(testInputPath, '64816L_amygdala_int.nii.gz')
  let testOutputFilePath = path.join(testOutputPath, '64816L_amygdala_int.nii.gz')

  const image = await readImageNode(testInputFilePath)
  const { outputImage } = await morphologicalContourInterpolationNode(image)
  await writeImageNode(outputImage, testOutputFilePath)

  testOutputFilePath = path.join(testOutputPath, '64816L_amygdala_int_axis_1.nii.gz')
  const { outputImage: outputImageAxis1 } = await morphologicalContourInterpolationNode(image, { axis: 1 })
  await writeImageNode(outputImageAxis1, testOutputFilePath)

  t.pass()
})
