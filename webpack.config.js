const MiniCssExtractPlugin = require('mini-css-extract-plugin');
const { CleanWebpackPlugin } = require('clean-webpack-plugin');
const UglifyJsPlugin = require('uglify-js-plugin');
const webpack = require("webpack");

const SCRIPTS = __dirname + "/webapp/";
const SCSS = __dirname + "/scss/";
const DEST = __dirname + "/docroot/scripts/";

module.exports = (env) => {

	const webpackConf = {
		entry: {
			'pdb-redo-style': SCSS + "pdb-redo-bootstrap.scss",
			'index': SCRIPTS + 'index.js'
		},

		module: {
			rules: [
				{
					test: /\.js$/,
					exclude: /node_modules/,
					use: {
						loader: "babel-loader",
						options: {
							presets: ['@babel/preset-env']
						}
					}
				},
				{
					test: /\.(eot|svg|ttf|woff(2)?)(\?v=\d+\.\d+\.\d+)?/,
					loader: 'file-loader',
					options: {
						name: '[name].[ext]',
						outputPath: 'fonts/',
						publicPath: '../fonts/'
					}
				},
				{
					test: /\.(png)/,
					loader: 'file-loader',
					options: {
						name: '[name].[ext]',
						outputPath: 'images/',
						publicPath: '../images/'
					}
				},
				{
					test: /\.s[ac]ss$/i,
					use: [
						MiniCssExtractPlugin.loader,
						'css-loader',
						'sass-loader'
					]
				}
			]
		},

		output: {
			path: DEST,
			filename: "./[name].js"
		},

		plugins: [
			new CleanWebpackPlugin({
				cleanOnceBeforeBuildPatterns: [
					'css/**/*',
					'scripts/**/*',
					'fonts/**/*',
					'!images*'
				]
			}),
			new webpack.ProvidePlugin({
				$: 'jquery',
				jQuery: 'jquery',
				CodeMirror: 'codemirror'
			}),
			new MiniCssExtractPlugin({
				filename: '../css/[name].css',
				chunkFilename: '../css/[id].css'
			})
		]
	};

    const PRODUCTIE = env != null && env.PRODUCTIE;

	if (PRODUCTIE) {
		webpackConf.mode = "production";
		webpackConf.plugins.push(new UglifyJsPlugin({parallel: 4}))
	} else {
		webpackConf.mode = "development";
		webpackConf.devtool = 'source-map';
		webpackConf.plugins.push(new webpack.optimize.AggressiveMergingPlugin())
	}

	return webpackConf;
};
